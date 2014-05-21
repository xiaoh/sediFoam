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
    labelList indexList,
    labelList cpuList
)
{
    label nprocs = Pstream::nProcs();

    labelList localIndexList(nprocs, 0);

    forAll(outList, listI)
    {
        outList[listI].setSize(indexList[listI]);
        forAll(outList[listI], i)
        {
            outList[listI][i] = 0*outList[listI][i];
        }
    }

    forAll(inList, i)
    {
        label processorI = cpuList[i];
        label& nIndex = localIndexList[processorI];
        outList[processorI][nIndex] = inList[i];
        nIndex ++;
    }
}

template <class DataType>
void softParticleCloud::flattenList
(
    List<DataType>& inList,
    DataType& outList
)
{
    label lengthList = 0.0;
    forAll(inList,listI)
    {
        lengthList += inList[listI].size();
    }

    outList.setSize(lengthList);

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

    labelList putLmpParticleNo(nprocs, 0);

    // calculate the number of particles in each LmpCpu
    vectorList putFoamDragList(size(), vector::zero);
    labelList putFoamTagList(size(), 0);
    labelList putFoamCpuIdList(size(), 0);

    labelList putLmpCpuIndexList(size(), 0);

    int i = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();

        putFoamTagList[i] = p.ptag();
        putFoamDragList[i] = FLocal[i];

        label putDragLmpCpuId = p.pLmpCpuId();
        putLmpCpuIndexList[i] = putDragLmpCpuId;
        putFoamCpuIdList[i] = myrank;
        putLmpParticleNo[putDragLmpCpuId] ++;
    }

    // set drag/tag list of particle in each LmpCpu
    // assemble them into lists of drag/tag/foamCpuId 
    // by the index of the processor
    List<vectorList> putDragListList(nprocs);
    List<labelList> putTagListList(nprocs);
    List<labelList> putFoamCpuIdListList(nprocs);
    assembleList<vectorList>
    (
        putFoamDragList,
        putDragListList,
        putLmpParticleNo,
        putLmpCpuIndexList
    );

    assembleList<labelList>
    (
        putFoamTagList,
        putTagListList,
        putLmpParticleNo,
        putLmpCpuIndexList
    );

    assembleList<labelList>
    (
        putFoamCpuIdList,
        putFoamCpuIdListList,
        putLmpParticleNo,
        putLmpCpuIndexList
    );

    // transpose the lists in each foamCpu to lmpCpu
    List<vectorList> putLmpDragListList(nprocs);
    List<labelList> putLmpTagListList(nprocs);
    List<labelList> putLmpFoamCpuIdListList(nprocs);

    transposeAmongProcs<vectorList> (putDragListList, putLmpDragListList);
    transposeAmongProcs<labelList> (putTagListList, putLmpTagListList);
    transposeAmongProcs<labelList>
    (
        putFoamCpuIdListList,
        putLmpFoamCpuIdListList
    );

    // flatten the lists obtained for each LmpCpu
    vectorList putLmpDragList;
    labelList putLmpTagList;
    labelList putLmpFoamCpuIdList;

    flattenList<vectorList> (putLmpDragListList, putLmpDragList);
    flattenList<labelList> (putLmpTagListList, putLmpTagList);
    flattenList<labelList> (putLmpFoamCpuIdListList, putLmpFoamCpuIdList);

    // transfer the data to the type that lammps can read
    int putListSize = putLmpDragList.size();
    double* putFLocalArray_ = new double [3*putListSize];
    int* putTagLocalArray_ = new int [putListSize];
    int* putFoamCpuIdLocalArray_ = new int [putListSize];

    for(label i = 0; i < putListSize; i++)
    {
        putFLocalArray_[3*i + 0] = putLmpDragList[i].x();
        putFLocalArray_[3*i + 1] = putLmpDragList[i].y();
        putFLocalArray_[3*i + 2] = putLmpDragList[i].z();
        putTagLocalArray_[i] = putLmpTagList[i];
        putFoamCpuIdLocalArray_[i] = putLmpFoamCpuIdList[i];
    }

    lammps_put_local_info
    (
        lmp_,
        putListSize,
        putFLocalArray_,
        putFoamCpuIdLocalArray_,
        putTagLocalArray_
    );

    delete [] putFLocalArray_;
    delete [] putTagLocalArray_;
    delete [] putFoamCpuIdLocalArray_;

    // Ask lammps to move certain steps forward
    lammps_step(lmp_, nstep);

    // Harvest the number of particles in each lmp cpu
    int lmpNLocal = lammps_get_local_n(lmp_);

    double* getXArrayLocal = new double [3*lmpNLocal];
    double* getVArrayLocal = new double [3*lmpNLocal];
    int* getFoamCpuIdArrayLocal = new int[lmpNLocal];
    int* getLmpCpuIdArrayLocal = new int[lmpNLocal];
    int* getTagArrayLocal = new int[lmpNLocal];

    lammps_get_local_info
    (
        lmp_,
        getXArrayLocal,
        getVArrayLocal,
        getFoamCpuIdArrayLocal,
        getLmpCpuIdArrayLocal,
        getTagArrayLocal
    );

    vectorList getXList(lmpNLocal, vector::zero);
    vectorList getVList(lmpNLocal, vector::zero);
    labelList getFoamCpuIdList(lmpNLocal, 0);
    labelList getLmpCpuIdList(lmpNLocal, 0);
    labelList getTagList(lmpNLocal, 0);

    labelList getFoamParticleNo(nprocs, 0);

    for(label i = 0; i < lmpNLocal; i++)
    {
        getXList[i] =
            vector
            (
                getXArrayLocal[3*i + 0],
                getXArrayLocal[3*i + 1],
                getXArrayLocal[3*i + 2]
            );

        getVList[i] =
            vector
            (
                getVArrayLocal[3*i + 0],
                getVArrayLocal[3*i + 1],
                getVArrayLocal[3*i + 2]
            );

        label foamCpuId = getFoamCpuIdArrayLocal[i];
        label lmpCpuId = getLmpCpuIdArrayLocal[i];
        getFoamCpuIdList[i] = foamCpuId;
        getLmpCpuIdList[i] = lmpCpuId;
        getTagList[i] = getTagArrayLocal[i];
        getFoamParticleNo[foamCpuId] ++;
    }

    delete [] getFoamCpuIdArrayLocal;
    delete [] getLmpCpuIdArrayLocal;
    delete [] getTagArrayLocal;
    delete [] getXArrayLocal;
    delete [] getVArrayLocal;

    List<vectorList> getXListList(nprocs);
    List<vectorList> getVListList(nprocs);
    List<labelList> getFoamCpuIdListList(nprocs);
    List<labelList> getLmpCpuIdListList(nprocs);
    List<labelList> getTagListList(nprocs);
    assembleList<vectorList>
    (
        getXList,
        getXListList,
        getFoamParticleNo,
        getFoamCpuIdList
    );

    assembleList<vectorList>
    (
        getVList,
        getVListList,
        getFoamParticleNo,
        getFoamCpuIdList
    );

    assembleList<labelList>
    (
        getFoamCpuIdList,
        getFoamCpuIdListList,
        getFoamParticleNo,
        getFoamCpuIdList
    );

    assembleList<labelList>
    (
        getLmpCpuIdList,
        getLmpCpuIdListList,
        getFoamParticleNo,
        getFoamCpuIdList
    );

    assembleList<labelList>
    (
        getTagList,
        getTagListList,
        getFoamParticleNo,
        getFoamCpuIdList
    );


    // transpose the lists in each LmpCpu to OfCpu
    List<vectorList> getFoamXListList(nprocs);
    List<vectorList> getFoamVListList(nprocs);
    List<labelList> getFoamFoamCpuIdListList(nprocs);
    List<labelList> getFoamLmpCpuIdListList(nprocs);
    List<labelList> getFoamTagListList(nprocs);

    transposeAmongProcs<vectorList> (getXListList, getFoamXListList);
    transposeAmongProcs<vectorList> (getVListList, getFoamVListList);
    transposeAmongProcs<labelList>
    (
        getFoamCpuIdListList,
        getFoamFoamCpuIdListList
    );
    transposeAmongProcs<labelList>
    (
        getLmpCpuIdListList,
        getFoamLmpCpuIdListList
    );
    transposeAmongProcs<labelList> (getTagListList, getFoamTagListList);

    vectorList getFoamXList;
    vectorList getFoamVList;
    labelList getFoamFoamCpuIdList;
    labelList getFoamLmpCpuIdList;
    labelList getFoamTagList;

    flattenList<vectorList> (getFoamXListList, getFoamXList);
    flattenList<vectorList> (getFoamVListList, getFoamVList);
    flattenList<labelList> (getFoamFoamCpuIdListList, getFoamFoamCpuIdList);
    flattenList<labelList> (getFoamLmpCpuIdListList, getFoamLmpCpuIdList);
    flattenList<labelList> (getFoamTagListList, getFoamTagList);

    int getListSize = getFoamXList.size();
    double* getXFoamArrayLocal = new double [3*getListSize];
    double* getVFoamArrayLocal = new double [3*getListSize];

    for(label i = 0; i < getListSize; i++)
    {
        label j = 0;
        for(j = 0; j < getListSize; j++)
        {
            if (putFoamTagList[i] == getFoamTagList[j])
            {
                break;
            }
        }

        getXFoamArrayLocal[3*i + 0] = getFoamXList[j].x();
        getXFoamArrayLocal[3*i + 1] = getFoamXList[j].y();
        getXFoamArrayLocal[3*i + 2] = getFoamXList[j].z();

        getVFoamArrayLocal[3*i + 0] = getFoamVList[j].x();
        getVFoamArrayLocal[3*i + 1] = getFoamVList[j].y();
        getVFoamArrayLocal[3*i + 2] = getFoamVList[j].z();
        
        lmpCpuIdLocal[i] = getFoamLmpCpuIdList[j];
    }

    // This a temperary work-around. Assemble array to vectors.
    assembleVectors(XLocal, getXFoamArrayLocal);
    assembleVectors(VLocal, getVFoamArrayLocal);

    delete [] getXFoamArrayLocal;
    delete [] getVFoamArrayLocal;

} // Job done; Proceed to next fluid calculation step.


void softParticleCloud::writeFields() const
{
    softParticle::writeFields(*this);
}

} // namespace Foam


// ************************************************************************* //
