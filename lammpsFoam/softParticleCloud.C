/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "softParticleCloud.H"
#include "processorPolyPatch.H"
#include "vectorList.H"
#include <string.h>
#include "mpi.h"
using std::string;

namespace Foam
{

    // defineParticleTypeNameAndDebug(softParticle, 0);
    defineTemplateTypeNameAndDebug(Cloud<softParticle>, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void softParticleCloud::initLammps()
{
    int i, N;

    string lamString;

    lamString =  string("in.lammps");
    N = lamString.length();

    char* lmpInputScript = new char[N + 1];
    lmpInputScript[N] = '\0';

    for (i = 0; i < N; i++)
        lmpInputScript[i] = lamString[i];

    MPI_Comm commLammps;
    // int iLammps = 1;
    //    MPI_Comm_split(MPI_COMM_WORLD,iLammps,0,);
    MPI_Comm_dup(MPI_COMM_WORLD, &commLammps);

    lmp_ = new LAMMPS(0,NULL,commLammps);

    FILE* fp = NULL;

    if (Pstream::master())
    {
        // Open in lammps input script
        fp = fopen(lmpInputScript,"r");
        if (fp == NULL)
        {
            printf("initLammps::ERROR: Could not open LAMMPS input script.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Scan through the file and run the file line by line
    Info<< "Reading Lammps inputfile (in.lammps) ..." << endl;

    lammps_sync(lmp_);

    int n = 0;
    char line[1024];

    while (1)
    {
        if (Pstream::master())
        {
            if (fgets(line,1024,fp) == NULL)
            {
                n = 0;
            }
            else
            {
                n = strlen(line) + 1;
            }

            if (n == 0) fclose(fp);
        }

        Pstream::scatter(n);
        if (n == 0) break;

        // The only raw MPI function in foam setcion.
        MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
        lmp_->input->one(line);

        string inputLine = line;
        if (inputLine.find("timestep",0) != string::npos)
        {
            Info<< "Timestep specified in Lammps input as follows: \n"
                << " --- " << line << " --- ";
            adjustLampTimestep();
        }
    }

    Info<< "Finished reading Lammps inputfile." << endl;
    // First, get no. of particles
    nGlobal_ = lammps_get_global_n(lmp_);
    Info<< "FOAM reported # of particles according to Lammps: "
        << nGlobal_ << endl;

    // Setup temporary space for holding x & v for sending to Lammps.
    // Note: global communication of particles is involved in the first step
    xArray_ = new double [3*nGlobal_];
    vArray_ = new double [3*nGlobal_];
    dArray_ = new double [nGlobal_];
    fArray_ = new double [3*nGlobal_];

    rhoArray_ = new double [nGlobal_];
    tagArray_ = new int [nGlobal_];
    lmpCpuIdArray_ = new int [nGlobal_];
    typeArray_ = new int [nGlobal_];

    // xArray_ etc. are local to Lammps processor
    lammps_get_initial_info
    (
        lmp_,
        xArray_,
        vArray_,
        dArray_,
        rhoArray_,
        tagArray_,
        lmpCpuIdArray_,
        typeArray_
    );

    initConstructParticles
    (
        nGlobal_,
        xArray_,
        vArray_,
        dArray_,
        rhoArray_,
        tagArray_,
        lmpCpuIdArray_,
        typeArray_
    );

    lammps_step(lmp_, 0);
}


void softParticleCloud::adjustLampTimestep()
{
    // Time step suggested in Lammps input file
    double dtLampIn = lammps_get_timestep(lmp_);

    // How many Lammps/solid steps (dtS) per fluid step (dtF)
    // Rounded upward (toward ceiling) to integer
    double dnSub = round(this->runTime().deltaT().value()/dtLampIn);
    if (dnSub == 0) dnSub++;

    solidStepsPerDt_ = scalar(dnSub);  // This is the real # of substeps.
    solidStepsPerDt_ =
      (int(solidStepsPerDt_)/int(subCycles_) + 1)*int(subCycles_);

    // Adjusted Lammps timestep
    scalar dtLampAdj = scalar(this->runTime().deltaT().value()/dnSub);

    Pstream::scatter(dtLampAdj);
    lammps_set_timestep(lmp_, dtLampAdj);

    if (subCycles_ >= solidStepsPerDt_)
    {
        subCycles_ = solidStepsPerDt_;
        nXtra_ = 0;
        subSteps_ = 1;
    }
    else
    {
        subSteps_ = int(solidStepsPerDt_)/int(subCycles_);
        nXtra_ = int(solidStepsPerDt_)%int(subCycles_);
        if (nXtra_ != 0)
        FatalErrorIn
        (
            "softParticleCloud::adjustLampTimestep() "
        )   << "Time step adjustment error. "
            << "Subcycles: " << subCycles_
            << "Steps per cycle: " <<  subSteps_
            << "Extra steps: " <<   nXtra_
            << abort(FatalError);
    }

    Pstream::scatter(subCycles_);
    Pstream::scatter(nXtra_);
    Pstream::scatter(subSteps_);

    Info<< "Lammps timestep adjusted from " << dtLampIn
        << " to " << dtLampAdj << endl;
    Info<< "Adjusted according to: " << int(solidStepsPerDt_)
        << " solid steps per fluid step." << endl;
    Info<< "Sub-cycle report: " << subCycles_ << " sub-cycles with "
        << subSteps_ << " steps per cycle and "
        << nXtra_ << " extra steps." << endl;
}


// Construct particles in FOAM from Lammps data:
void softParticleCloud::initConstructParticles
(
   int nGlobal_,
   double* x,
   double* v,
   double* d,
   double* rho,
   int* tag,
   int* lmpCpuId,
   int* type
)
{
    int offset;
    for (int i = 0; i < nGlobal_; i++)
    {
        offset = 3*i;
        vector pos = vector
        (
            x[offset + 0],
            x[offset + 1],
            x[offset + 2]
        );

        vector velo = vector
        (
            v[offset + 0],
            v[offset + 1],
            v[offset + 2]
        );

        label cellI = mesh_.findCell(pos);

        if (cellI >= 0) ++nLocal_;

        label gCellI = cellI;

        reduce(gCellI, maxOp<label>());

        if (gCellI < 0)
        {
            FatalErrorIn
            (
                "softParticleCloud::initConstructParticles "
            )   << "Particle leaving domain. \n "
                << "Proc #: " << Pstream::myProcNo() << endl
                << "Particle #: " << i << endl
                << "Projected particle position: " << pos << endl
                << abort(FatalError);
        }

        scalar ds = scalar(d[i]);
        scalar rhos = scalar(rho[i]);
        label tags = int(tag[i]);
        label lmpCpuIds = int(lmpCpuId[i]);
        label types = int(type[i]);

        // create a new softParticle when it is in the current processor
        // but the computer is running much slower than before.
        if (cellI >= 0)
        {
            softParticle* ptr =
                new softParticle
                (
                    pMesh(),
                    pos,
                    cellI,
                    ds,
                    velo,
                    rhos,
                    tags,
                    lmpCpuIds,
                    types
                );

            if (debug)
            {
                Pout<< "position is:" << pos << endl;
                Pout<< "cell is:" << cellI << endl;
                Pout<< "foam tag is:" << tags << endl;
                Pout<< "lammps CPU id is:" << lmpCpuId << endl;
                Pout<< "type is:" << types << endl;
            }
            addParticle(ptr);
        }
    }
    Pout<< " After initialization, I have " << nLocal_
        << " particles locally." << endl;
}


// Wrap up Lammps, i.e. delete the pointer.
void softParticleCloud::finishLammps()
{
    if (Pstream::master())
        delete lmp_;

    delete [] xArray_;
    delete [] vArray_;
    delete [] dArray_;
    delete [] fArray_;
    delete [] rhoArray_;
    delete [] tagArray_;
    delete [] lmpCpuIdArray_;
    delete [] typeArray_;
}


// The two function will be eliminated in the future for efficiency
// reasons. This is a work around.

// This is not necessary as vector in OF is contiguous!
// Will remove later.
void softParticleCloud::flattenVectors
(
    const vector* vec,
    double* array
)
{
    for(label i = 0; i < nGlobal_; i++)
    {
        array[3*i + 0] = vec[i].x();
        array[3*i + 1] = vec[i].y();
        array[3*i + 2] = vec[i].z();
    }
}

void softParticleCloud::flattenVectors
(
    const vectorList vec,
    double* array
)
{
    for(label i = 0; i < nGlobal_; i++)
    {
        array[3*i + 0] = vec[i].x();
        array[3*i + 1] = vec[i].y();
        array[3*i + 2] = vec[i].z();

        if (debug)
        {
            Pout<< "array [" << 3*i + 0 << "] is: "
                << array[3*i + 0] << endl;
            Pout<< "array [" << 3*i + 1 << "] is: "
                << array[3*i + 1] << endl;
            Pout<< "array [" << 3*i + 2 << "] is: "
                << array[3*i + 2] << endl;
        }
    }
}


// Reverse operation of flattenVectors;
void softParticleCloud::assembleVectors
(
    vector* vec,
    const double* array
)
{
    for(label i = 0; i < size(); i++)
    {
        vec[i] = vector
        (
            array[3*i + 0],
            array[3*i + 1],
            array[3*i + 2]
        );
        if (debug)
        {
            Pout<< "vec [" << i << "] is: " << vec[i] << endl;
        }
    }
}


template <class DataType>
void softParticleCloud::assembleList
(
    DataType& inList,
    List<DataType>& outList,
    labelList listSizeList,
    labelList cpuIdList
)
{
    label nprocs = Pstream::nProcs();

    labelList localIndexList(nprocs, 0);

    forAll(outList, listI)
    {
        outList[listI].setSize(listSizeList[listI]);
    }

    forAll(inList, i)
    {
        label processorI = cpuIdList[i];
        label& nIndex = localIndexList[processorI];
        outList[processorI][nIndex] = inList[i];
        nIndex++;
    }
}

template <class DataType>
void softParticleCloud::flattenList
(
    List<DataType>& inList,
    DataType& outList
)
{
    label listLength = 0.0;
    forAll(inList,listI)
    {
        listLength += inList[listI].size();
    }

    outList.setSize(listLength);

    int i = 0;
    forAll(inList,listI)
    {
        forAll(inList[listI],memberI)
        {
            outList[i] = inList[listI][memberI];
            i++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
softParticleCloud::softParticleCloud
(
    const volPointInterpolation& vpi,
    const volVectorField& U,
    const volScalarField& p,
    volVectorField& Ue,
    dimensionedScalar nu,
    volScalarField& alpha,
    IOdictionary& cloudProperties
)
:
    Cloud<softParticle>(U.mesh(), "softParticleCloud"),
    nGlobal_(0),
    runTime_(U.time()),
    mesh_(U.mesh()),
    nu_(nu),
    volPointInterpolation_(vpi),
    cloudProperties_(cloudProperties),
    Ue_(Ue),
    U_(U),
    pf_(p),
    gamma_(alpha)
{
    subCycles_ = readScalar(cloudProperties_.lookup("subCycles"));

    nLocal_ = 0;

    // This is particularly for continuation run. After reading in
    // diameters, refresh gamma field.
    // ===> If continuation run, ask Lammps to read restart file.
    // To be implemented.

    initLammps();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

softParticleCloud::~softParticleCloud()
{
    finishLammps();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Change the positions and velocities of all the particles in this
// cloud according to information provided by Lammps.
void softParticleCloud::setPositionVeloCpuId
(
    vector XLocal[],
    vector VLocal[],
    int lmpCpuIdLocal[]
)
{
    int i = 0;
    nLocal_ = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();

        // Update position:
        p.positionOld() = p.position();
        p.moveU() = (XLocal[i] - p.position())/mesh_.time().deltaTValue();

        // Update velocity:
        p.U() = VLocal[i];

        p.pLmpCpuId() = lmpCpuIdLocal[i];
        nLocal_++;
    }
    Pout<< " After movement, I have " << nLocal_
        << " particles locally." << endl;
}

//- Set the particle cell index after the particle
//  moves across the processor boundary
void softParticleCloud::setPositionCell()
{
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter
    )
    {
        softParticle& p = pIter();

        // Update cell number:
        p.cell() = mesh_.findCell(p.position());
    }
}

//  Temp placement;

// transpose information on processors

template<class DataType>
void softParticleCloud:: transposeAmongProcs
(
    List<DataType>& toEveryone,
    List<DataType>& fromEveryone
)
{

    label nprocs = Pstream::nProcs();
    label myrank = Pstream::myProcNo();

    if (toEveryone.size() != nprocs)
    {
        FatalErrorIn
	    (
	     "softParticleCloud::transAmongProcs() "
	    )   << "List size not equal to no. of procs \n "
            << "Proc #: " << myrank
            << abort(FatalError);
    }

    PstreamBuffers pBufs(Pstream::nonBlocking);

    // send to everyone else with UOPstream (output)
    for(int i = 0; i < nprocs; i++)
    {
        if (i != myrank)
        {
            UOPstream toEveryoneStream(i, pBufs);
            toEveryoneStream << toEveryone[i];
        }
    }

    // finish sending and block threads here
    pBufs.finishedSends();

    // receive from everyone else; copy local data residing on my proc
    for(int i = 0; i < nprocs; i++)
    {
        if (i != myrank)
        {
            UIPstream fromEveryoneStream(i, pBufs);
            fromEveryoneStream >> fromEveryone[i];
        }
        else
        {
            fromEveryone[i] = toEveryone[i];
        }
    }
}

// #define DEBUG_EVOLVE

// Call Lammps and evolve forward some steps.  Given particle
// positions, velocities, and drag.
// This is a wrapper for the Lammps functions which actually does
// the job.
void  softParticleCloud::lammpsEvolveForward
(
    vector* XLocal,
    vector* VLocal,
    int* lmpCpuIdLocal,
    vectorList FLocal,
    int nstep
)
{
    label nprocs = Pstream::nProcs();
    label myrank = Pstream::myProcNo();

    // Start putting information to LAMMPS
    labelList lmpParticleNo(nprocs, 0);

    // Calculate the number of particles in each LmpCpu
    // and obtain the drag force
    vectorList fromFoamDragList(size(), vector::zero);
    labelList fromFoamFoamCpuIdList(size(), myrank);
    labelList fromFoamTagList(size(), 0);

    labelList assembleLmpCpuIdList(size(), 0);

    int i = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();

        fromFoamDragList[i] = FLocal[i];
        fromFoamTagList[i] = p.ptag();

        label lmpCpuId = p.pLmpCpuId();
        assembleLmpCpuIdList[i] = lmpCpuId;
        lmpParticleNo[lmpCpuId] ++;
    }

    // Set drag/tag list of particle in each LmpCpu
    // assemble them into lists of drag/tag/foamCpuId 
    // by the index of the processor
    List<vectorList> fromFoamDragListList(nprocs);
    List<labelList> fromFoamFoamCpuIdListList(nprocs);
    List<labelList> fromFoamTagListList(nprocs);
    assembleList<vectorList>
    (
        fromFoamDragList,
        fromFoamDragListList,
        lmpParticleNo,
        assembleLmpCpuIdList
    );

    assembleList<labelList>
    (
        fromFoamFoamCpuIdList,
        fromFoamFoamCpuIdListList,
        lmpParticleNo,
        assembleLmpCpuIdList
    );

    assembleList<labelList>
    (
        fromFoamTagList,
        fromFoamTagListList,
        lmpParticleNo,
        assembleLmpCpuIdList
    );

    // transpose the lists in each foamCpu to lmpCpu
    List<vectorList> toLmpDragListList(nprocs);
    List<labelList> toLmpFoamCpuIdListList(nprocs);
    List<labelList> toLmpTagListList(nprocs);

    transposeAmongProcs<vectorList> (fromFoamDragListList, toLmpDragListList);
    transposeAmongProcs<labelList>
    (
        fromFoamFoamCpuIdListList,
        toLmpFoamCpuIdListList
    );
    transposeAmongProcs<labelList> (fromFoamTagListList, toLmpTagListList);

    // flatten the lists obtained for each LmpCpu
    vectorList toLmpDragList;
    labelList toLmpFoamCpuIdList;
    labelList toLmpTagList;

    flattenList<vectorList> (toLmpDragListList, toLmpDragList);
    flattenList<labelList> (toLmpFoamCpuIdListList, toLmpFoamCpuIdList);
    flattenList<labelList> (toLmpTagListList, toLmpTagList);

    // transfer the data to the type that lammps can read
    int toLmpListSize = toLmpDragList.size();
    double* toLmpFLocalArray_ = new double [3*toLmpListSize];
    int* toLmpFoamCpuIdLocalArray_ = new int [toLmpListSize];
    int* toLmpTagLocalArray_ = new int [toLmpListSize];

    for(label i = 0; i < toLmpListSize; i++)
    {
        toLmpFLocalArray_[3*i + 0] = toLmpDragList[i].x();
        toLmpFLocalArray_[3*i + 1] = toLmpDragList[i].y();
        toLmpFLocalArray_[3*i + 2] = toLmpDragList[i].z();
        toLmpFoamCpuIdLocalArray_[i] = toLmpFoamCpuIdList[i];
        toLmpTagLocalArray_[i] = toLmpTagList[i];
    }

    lammps_put_local_info
    (
        lmp_,
        toLmpListSize,
        toLmpFLocalArray_,
        toLmpFoamCpuIdLocalArray_,
        toLmpTagLocalArray_
    );

    delete [] toLmpFLocalArray_;
    delete [] toLmpFoamCpuIdLocalArray_;
    delete [] toLmpTagLocalArray_;

    // Ask lammps to move certain steps forward
    lammps_step(lmp_, nstep);

    // Start getting information from LAMMPS
    // Harvest the number of particles in each lmp cpu
    int lmpNLocal = lammps_get_local_n(lmp_);

    // Harvest more infomation from each lmp cpu
    double* fromLmpXArrayLocal = new double [3*lmpNLocal];
    double* fromLmpVArrayLocal = new double [3*lmpNLocal];
    int* fromLmpFoamCpuIdArrayLocal = new int[lmpNLocal];
    int* fromLmpLmpCpuIdArrayLocal = new int[lmpNLocal];
    int* fromLmpTagArrayLocal = new int[lmpNLocal];

    lammps_get_local_info
    (
        lmp_,
        fromLmpXArrayLocal,
        fromLmpVArrayLocal,
        fromLmpFoamCpuIdArrayLocal,
        fromLmpLmpCpuIdArrayLocal,
        fromLmpTagArrayLocal
    );

    // Transform the data obtained from lammps to openfoam format
    vectorList fromLmpXList(lmpNLocal, vector::zero);
    vectorList fromLmpVList(lmpNLocal, vector::zero);
    labelList fromLmpFoamCpuIdList(lmpNLocal, 0);
    labelList fromLmpLmpCpuIdList(lmpNLocal, myrank);
    labelList fromLmpTagList(lmpNLocal, 0);

    labelList foamParticleNo(nprocs, 0);

    for(label i = 0; i < lmpNLocal; i++)
    {
        fromLmpXList[i] =
            vector
            (
                fromLmpXArrayLocal[3*i + 0],
                fromLmpXArrayLocal[3*i + 1],
                fromLmpXArrayLocal[3*i + 2]
            );

        fromLmpVList[i] =
            vector
            (
                fromLmpVArrayLocal[3*i + 0],
                fromLmpVArrayLocal[3*i + 1],
                fromLmpVArrayLocal[3*i + 2]
            );

        label foamCpuId = fromLmpFoamCpuIdArrayLocal[i];
        fromLmpFoamCpuIdList[i] = foamCpuId;
        fromLmpTagList[i] = fromLmpTagArrayLocal[i];
        foamParticleNo[foamCpuId] ++;
    }

    delete [] fromLmpXArrayLocal;
    delete [] fromLmpVArrayLocal;
    delete [] fromLmpFoamCpuIdArrayLocal;
    delete [] fromLmpLmpCpuIdArrayLocal;
    delete [] fromLmpTagArrayLocal;

    // Separate the information from LAMMPS into  lists of 
    // different processors
    List<vectorList> fromLmpXListList(nprocs);
    List<vectorList> fromLmpVListList(nprocs);
    List<labelList> fromLmpFoamCpuIdListList(nprocs);
    List<labelList> fromLmpLmpCpuIdListList(nprocs);
    List<labelList> fromLmpTagListList(nprocs);
    assembleList<vectorList>
    (
        fromLmpXList,
        fromLmpXListList,
        foamParticleNo,
        fromLmpFoamCpuIdList
    );

    assembleList<vectorList>
    (
        fromLmpVList,
        fromLmpVListList,
        foamParticleNo,
        fromLmpFoamCpuIdList
    );

    assembleList<labelList>
    (
        fromLmpFoamCpuIdList,
        fromLmpFoamCpuIdListList,
        foamParticleNo,
        fromLmpFoamCpuIdList
    );

    assembleList<labelList>
    (
        fromLmpLmpCpuIdList,
        fromLmpLmpCpuIdListList,
        foamParticleNo,
        fromLmpFoamCpuIdList
    );

    assembleList<labelList>
    (
        fromLmpTagList,
        fromLmpTagListList,
        foamParticleNo,
        fromLmpFoamCpuIdList
    );

    // Transpose the lists in each LmpCpu to FoamCpu
    List<vectorList> toFoamXListList(nprocs);
    List<vectorList> toFoamVListList(nprocs);
    List<labelList> toFoamFoamCpuIdListList(nprocs);
    List<labelList> toFoamLmpCpuIdListList(nprocs);
    List<labelList> toFoamTagListList(nprocs);

    transposeAmongProcs<vectorList> (fromLmpXListList, toFoamXListList);
    transposeAmongProcs<vectorList> (fromLmpVListList, toFoamVListList);
    transposeAmongProcs<labelList>
    (
        fromLmpFoamCpuIdListList,
        toFoamFoamCpuIdListList
    );
    transposeAmongProcs<labelList>
    (
        fromLmpLmpCpuIdListList,
        toFoamLmpCpuIdListList
    );
    transposeAmongProcs<labelList> (fromLmpTagListList, toFoamTagListList);

    // Collect all the information transfered from other processors in lists
    // and combine them into one list
    vectorList toFoamXList;
    vectorList toFoamVList;
    labelList toFoamFoamCpuIdList;
    labelList toFoamLmpCpuIdList;
    labelList toFoamTagList;

    flattenList<vectorList> (toFoamXListList, toFoamXList);
    flattenList<vectorList> (toFoamVListList, toFoamVList);
    flattenList<labelList> (toFoamFoamCpuIdListList, toFoamFoamCpuIdList);
    flattenList<labelList> (toFoamLmpCpuIdListList, toFoamLmpCpuIdList);
    flattenList<labelList> (toFoamTagListList, toFoamTagList);

    // Assign the position & velocity & lmpCpuId to the particle in OpenFOAM
    int toListSize = toFoamXList.size();

    for(label i = 0; i < toListSize; i++)
    {
        label j = 0;
        for(j = 0; j < toListSize; j++)
        {
            if (fromFoamTagList[i] == toFoamTagList[j])
            {
                break;
            }
        }

        // position
        XLocal[i] = 
            vector
            (
                toFoamXList[j].x(),
                toFoamXList[j].y(),
                toFoamXList[j].z()
            );

        // velocity
        VLocal[i] = 
            vector
            (
                toFoamVList[j].x(),
                toFoamVList[j].y(),
                toFoamVList[j].z()
            );
        
        // lmpCpuId
        lmpCpuIdLocal[i] = toFoamLmpCpuIdList[j];
    }

} // Job done; Proceed to next fluid calculation step.


void softParticleCloud::writeFields() const
{
    softParticle::writeFields(*this);
}

} // namespace Foam


// ************************************************************************* //
