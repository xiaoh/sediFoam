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
#include "tensorList.H"
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

    label nprocs = Pstream::nProcs();
    label myrank = Pstream::myProcNo();

    int* npArray = new int [nprocs];

    lammps_get_initial_np(lmp_, npArray);

    Info<< "creating new arrays..." << endl;
    Info<< "execution time is: " << runTime_.elapsedCpuTime() << endl;

    Pout<< "npArray is: " << npArray[myrank] << endl;
    int nLocal = npArray[myrank];

    xArray_ = new double [3*nLocal];
    vArray_ = new double [3*nLocal];
    dArray_ = new double [nLocal];
    fArray_ = new double [3*nLocal];

    rhoArray_ = new double [nLocal];
    tagArray_ = new int [nLocal];
    lmpCpuIdArray_ = new int [nLocal];
    typeArray_ = new int [nLocal];

    Info<< "getting initial info of particles..." << endl;
    Info<< "execution time is: " << runTime_.elapsedCpuTime() << endl;

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

    Info<< "constructing particles..." << endl;
    Info<< "execution time is: " << runTime_.elapsedCpuTime() << endl;

    initConstructParticles
    (
        nLocal,
        xArray_,
        vArray_,
        dArray_,
        rhoArray_,
        tagArray_,
        lmpCpuIdArray_,
        typeArray_
    );


    Info<< "before moving..." << endl;
    Info<< "execution time is: " << runTime_.elapsedCpuTime() << endl;

    softParticle::trackingData td0(*this);
    Cloud<softParticle>::move(td0, mesh_.time().deltaTValue());

    Info<< "execution time is: " << runTime_.elapsedCpuTime() << endl;

    lammps_step(lmp_, 0);

    double* lmpLocalBox = new double [6];
    lammps_get_local_domain(lmp_, lmpLocalBox);

    for (int i = 0; i < 6; i++)
    {
        lmpLocalBoxList_[myrank].component(i) = lmpLocalBox[i];
    }

    for(int i = 0; i < nprocs; i++)
    {
        reduce(lmpLocalBoxList_[i], sumOp<tensor>());
    }

    delete [] npArray;
    delete [] lmpLocalBox;
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
      (int(solidStepsPerDt_)/int(subCycles_))*int(subCycles_);

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
   int nLocal,
   double* x,
   double* v,
   double* d,
   double* rho,
   int* tag,
   int* lmpCpuId,
   int* type
)
{
    // if (Pstream::master())
    {
        int offset;
        for (int i = 0; i < nLocal; i++)
        {
            offset = 3*i;
            vector pos = mesh_.C()[0];
            label cellI = 0;

            vector velo = vector
            (
                v[offset + 0],
                v[offset + 1],
                v[offset + 2]
            );

            scalar ds = scalar(d[i]);
            scalar rhos = scalar(rho[i]);
            label tags = int(tag[i]);
            label lmpCpuIds = int(lmpCpuId[i]);
            label types = int(type[i]);

            // create a new softParticle when it is in the current processor
            // but the computer is running much slower than before.
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

        label i = 0;
        for
        (
            softParticleCloud::iterator pIter = begin();
            pIter != end();
            ++pIter, ++i
        )
        {
            offset = 3*i;

            vector pos = vector
            (
                x[offset + 0],
                x[offset + 1],
                x[offset + 2]
            );

            softParticle& p = pIter();

            // Update position:
            p.positionOld() = p.position();
            p.moveU() = (pos - p.position())/mesh_.time().deltaTValue();
        }
    }
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
    cloudProperties_(cloudProperties),
    Ue_(Ue),
    U_(U),
    pf_(p),
    gamma_(alpha),
    cpuTimeSplit_(6, 0.0)
{
    label nprocs = Pstream::nProcs();

    subCycles_ = readScalar(cloudProperties_.lookup("subCycles"));

    transposeNbrOnly_ =
        cloudProperties_.lookupOrDefault<Switch>("transposeNbrOnly", false);

    // Initialize the setup of adding and deleting particles
    addParticleOption_ = cloudProperties_.lookupOrDefault("addParticle", 0);
    deleteParticleOption_ = cloudProperties_.lookupOrDefault("deleteParticle", 0);
    deleteBeforeAddFlag_ = cloudProperties_.lookupOrDefault("deleteBeforeAdd", 0);

    gravity_ = cloudProperties_.lookupOrDefault("g", vector(0,-9.8,0));


    // initial the setup when adding particle
    if (addParticleOption_ > 0)
    {
        addParticleTimeStep_ =
            readScalar(cloudProperties_.lookup("addParticleTimeStep"));
        randomPerturb_ =
            readScalar(cloudProperties_.lookup("randomPerturb"));
        addParticleBox_ =
            cloudProperties_.lookupOrDefault("addParticleBox", tensor::zero);
        addParticleBoxEccentricity_ =
            cloudProperties_.lookupOrDefault("eccentricity", vector::zero);
        addParticleInfo_ = cloudProperties_.lookup("addParticleInfo");
        inletBox_ = cloudProperties_.lookupOrDefault("inletBox", tensor::zero);
        timeToAddParticle_ = addParticleTimeStep_;
        totalAdd_ = 0;
    }
    else
    {
        addParticleTimeStep_ = 1/SMALL;
        randomPerturb_ = 0.0;
        addParticleBox_ = tensor::zero;
        addParticleInfo_ = vector::zero;
        inletBox_ = tensor::zero;
        addParticleBoxEccentricity_ = vector::zero;
        timeToAddParticle_ = addParticleTimeStep_;
        totalAdd_ = 0;
    }

    // initial the setup when deleting particle
    if (deleteParticleOption_ > 0)
    {
        deleteParticleBox_ =
            cloudProperties_.lookupOrDefault("deleteParticleBox", tensor::zero);
        totalDelete_ = 0;
    }
    else
    {
        deleteParticleBox_ = tensor::zero;
        totalDelete_ = 0;
    }

    // initial the setup when clear the region where the particles are added
    if (deleteBeforeAddFlag_ > 0)
    {
        clearInitialBox_ =
            cloudProperties_.lookupOrDefault("clearInitialBox", tensor::zero);
        totalDeleteBeforeAdd_ = 0;
    }
    else
    {
        clearInitialBox_ = tensor::zero;
        totalDeleteBeforeAdd_ = 0;
    }

    nLocal_ = 0;

    lmpLocalBoxList_.setSize(nprocs);
    forAll(lmpLocalBoxList_, i)
    {
        lmpLocalBoxList_[i] = tensor::zero;
    }

    findAddParticleCells();

    // This is particularly for continuation run. After reading in
    // diameters, refresh gamma field.
    // ===> If continuation run, ask Lammps to read restart file.
    // To be implemented.

    initLammps();
    Info<< "initialization finished!" << endl;
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
        p.UOld() = p.U();
        p.U() = VLocal[i];

        p.pLmpCpuId() = lmpCpuIdLocal[i];
        nLocal_++;
    }
    // Pout<< " After movement, I have " << nLocal_
    //     << " particles locally." << endl;
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

    if (transposeNbrOnly_)
    {
        PstreamBuffers pBufsNbr(Pstream::nonBlocking);
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UOPstream toNbr(procPatch.neighbProcNo(), pBufsNbr);
                toNbr << toEveryone[procPatch.neighbProcNo()];
            }
        }

        pBufsNbr.finishedSends();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UIPstream fromNbr(procPatch.neighbProcNo(), pBufsNbr);
                fromNbr >> fromEveryone[procPatch.neighbProcNo()];
            }
        }

        fromEveryone[myrank] = toEveryone[myrank];
    }
    else
    {
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
            else if (i == myrank)
            {
                fromEveryone[i] = toEveryone[i];
            }
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

    label nList = size();

    // Start putting information to LAMMPS
    labelList lmpParticleNo(nprocs, 0);

    scalar t0 = runTime_.elapsedCpuTime();
    // Calculate the number of particles in each LmpCpu
    // and obtain the drag force
    vectorList fromFoamDragList(nList, vector::zero);
    labelList fromFoamFoamCpuIdList(nList, myrank);
    labelList fromFoamTagList(nList, 0);

    labelList assembleLmpCpuIdList(nList, 0);

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

    cpuTimeSplit_[0] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

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

    cpuTimeSplit_[1] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

    // flatten the lists obtained for each LmpCpu
    label toLmpListSize = 0;
    forAll(toLmpTagListList, listI)
    {
        toLmpListSize += toLmpTagListList[listI].size();
    }
    vectorList toLmpDragList(toLmpListSize, vector::zero);
    labelList toLmpFoamCpuIdList(toLmpListSize, 0);
    labelList toLmpTagList(toLmpListSize, 0);

    flattenList<vectorList> (toLmpDragListList, toLmpDragList);
    flattenList<labelList> (toLmpFoamCpuIdListList, toLmpFoamCpuIdList);
    flattenList<labelList> (toLmpTagListList, toLmpTagList);

    cpuTimeSplit_[2] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

    addAndDeleteParticle();

    if (contiguous<vector>())
    {
        // Pout<< "contiguous vector is true." << endl;
        double* toLmpFLocalArray_ =
            reinterpret_cast <double *> (&(toLmpDragList.first()));

        int* toLmpFoamCpuIdLocalArray_ =
            reinterpret_cast <int *> (&(toLmpFoamCpuIdList.first()));

        int* toLmpTagLocalArray_ =
            reinterpret_cast <int *> (&(toLmpTagList.first()));

        lammps_put_local_info
        (
            lmp_,
            toLmpListSize,
            toLmpFLocalArray_,
            toLmpFoamCpuIdLocalArray_,
            toLmpTagLocalArray_
        );
    }
    else
    {
        double* toLmpFLocalArray_ = new double [3*toLmpListSize];
        int* toLmpFoamCpuIdLocalArray_ = new int [toLmpListSize];
        int* toLmpTagLocalArray_ = new int [toLmpListSize];

        // Pout<< "contiguous vector is false." << endl;
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
    }

    // delete [] toLmpFLocalArray_;
    // delete [] toLmpFoamCpuIdLocalArray_;
    // delete [] toLmpTagLocalArray_;

    cpuTimeSplit_[3] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

    // Ask lammps to move certain steps forward
    Info<< "LAMMPS evolving.. " << endl;
    lammps_step(lmp_, nstep);

    Info<< "finished moving the particles in LAMMPS." << endl;
    cpuTimeSplit_[4] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();
    // Start getting information from LAMMPS
    // Harvest the number of particles in each lmp cpu
    int lmpNLocal = lammps_get_local_n(lmp_);

    label lmpNGlobal = lmpNLocal;
    reduce(lmpNGlobal, sumOp<label>());

    Info<< "the number of particles in LAMMPS now is: " << lmpNGlobal << endl;

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

    Info<< "data transformed into each processor!" << endl;

    delete [] fromLmpXArrayLocal;
    delete [] fromLmpVArrayLocal;
    delete [] fromLmpFoamCpuIdArrayLocal;
    delete [] fromLmpLmpCpuIdArrayLocal;
    delete [] fromLmpTagArrayLocal;

    cpuTimeSplit_[4] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

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

    // Separate the information from LAMMPS into  lists of
    cpuTimeSplit_[0] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

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

    cpuTimeSplit_[1] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

    // Collect all the information transfered from other processors in lists
    // and combine them into one list
    vectorList toFoamXList(nList, vector::zero);
    vectorList toFoamVList(nList, vector::zero);
    labelList toFoamFoamCpuIdList(nList, 0);
    labelList toFoamLmpCpuIdList(nList, 0);
    labelList toFoamTagList(nList, 0);

    flattenList<vectorList> (toFoamXListList, toFoamXList);
    flattenList<vectorList> (toFoamVListList, toFoamVList);
    flattenList<labelList> (toFoamFoamCpuIdListList, toFoamFoamCpuIdList);
    flattenList<labelList> (toFoamLmpCpuIdListList, toFoamLmpCpuIdList);
    flattenList<labelList> (toFoamTagListList, toFoamTagList);

    // Assign the position & velocity & lmpCpuId to the particle in OpenFOAM
    labelList sortedFromFoamTag(nList,0);
    labelList sortedToFoamTag(nList,0);

    sortedOrder(fromFoamTagList, sortedFromFoamTag);
    sortedOrder(toFoamTagList, sortedToFoamTag);

    for(label i = 0; i < nList; i++)
    {
        label fromI = sortedFromFoamTag[i];
        label toI = sortedToFoamTag[i];

        // position
        XLocal[fromI] =
            vector
            (
                toFoamXList[toI].x(),
                toFoamXList[toI].y(),
                toFoamXList[toI].z()
            );

        // velocity
        VLocal[fromI] =
            vector
            (
                toFoamVList[toI].x(),
                toFoamVList[toI].y(),
                toFoamVList[toI].z()
            );

        // lmpCpuId
        lmpCpuIdLocal[fromI] = toFoamLmpCpuIdList[toI];
    }

    cpuTimeSplit_[2] += runTime_.elapsedCpuTime() - t0;
    t0 = runTime_.elapsedCpuTime();

    Info<< "LAMMPS evolving finished! .. " << endl;
} // Job done; Proceed to next fluid calculation step.


// Add new particles in OpenFOAM
void softParticleCloud::addNewParticles()
{
    Info << "Adding new particle... " << endl;

    label nprocs = Pstream::nProcs();
    label myrank = Pstream::myProcNo();

    int npOF = addParticleCellID_.size();

    labelList addedParticleNo(nprocs, 0);
    labelList assembleLmpCpuIdList(npOF);

    vectorList fromFoamAddPositionList(npOF);
    List<vectorList> fromFoamAddPositionListList(nprocs);
    labelList fromFoamAddTagList(npOF);
    List<labelList> fromFoamAddTagListList(nprocs);

    forAll(addParticleCellID_, i)
    {
        label cellI = addParticleCellID_[i];
        vector pos = mesh_.C()[cellI];
        fromFoamAddPositionList[i] = pos;
        fromFoamAddTagList[i] = maxTag_ + 1 + i;
        forAll(lmpLocalBoxList_, boxI)
        {
            if (pointInBox(pos, lmpLocalBoxList_[boxI]))
            {
                addedParticleNo[boxI] ++;
                assembleLmpCpuIdList[i] = boxI;
            }
        }
    }

    assembleList<vectorList>
    (
        fromFoamAddPositionList,
        fromFoamAddPositionListList,
        addedParticleNo,
        assembleLmpCpuIdList
    );

    assembleList<labelList>
    (
        fromFoamAddTagList,
        fromFoamAddTagListList,
        addedParticleNo,
        assembleLmpCpuIdList
    );

    List<vectorList> toLmpAddPositionListList(nprocs);
    List<labelList> toLmpAddTagListList(nprocs);

    transposeAmongProcs<vectorList> (fromFoamAddPositionListList, toLmpAddPositionListList);
    transposeAmongProcs<labelList> (fromFoamAddTagListList, toLmpAddTagListList);

    label toLmpListSize = 0;
    forAll(toLmpAddPositionListList, listI)
    {
        toLmpListSize += toLmpAddPositionListList[listI].size();
    }
    vectorList toLmpAddPositionList(toLmpListSize, vector::zero);
    labelList toLmpAddTagList(toLmpListSize, 0);

    flattenList<vectorList> (toLmpAddPositionListList, toLmpAddPositionList);
    flattenList<labelList> (toLmpAddTagListList, toLmpAddTagList);

    int npAdd = toLmpListSize;

    double* posArray = new double [3*npAdd];
    double* tagArray = new double [npAdd];

    double ds = addParticleInfo_[0];
    double rhos = addParticleInfo_[1];
    int types = int(addParticleInfo_[2]);

    Random perturbation(size());

    forAll(toLmpAddPositionList, i)
    {
        posArray[0+3*i] =
            toLmpAddPositionList[i].x()
          + randomPerturb_*(0.5 - perturbation.scalar01());
        posArray[1+3*i] =
            toLmpAddPositionList[i].y()
          + randomPerturb_*(0.5 - perturbation.scalar01());
        posArray[2+3*i] =
            toLmpAddPositionList[i].z()
          + randomPerturb_*(0.5 - perturbation.scalar01());
        tagArray[i] = toLmpAddTagList[i];
    }

    label npAddGlobal = npAdd;
    reduce(npAddGlobal, sumOp<label>());
    totalAdd_ += npAddGlobal;

    lammps_create_particle(lmp_, npAdd, posArray, tagArray, ds, rhos, types);
    delete [] posArray;
    delete [] tagArray;
}


//- Adding and deleting particles
void softParticleCloud::addAndDeleteParticle()
{
    // Try adding new particles
    if (addParticleOption_ > 0 && timeToAddParticle_ <= 0)
    {
        // Deleting particles before add
        if (deleteBeforeAddFlag_ == 1)
        {
            int size = deleteBeforeAddList_.size();
            int* deleteList = new int [size];
            int nDelete = 0;
            for (int i = 0; i < size; i++)
            {
                // TODO: this may have problem for parallel computing
                deleteList[i] = deleteBeforeAddList_[i];
                if (deleteList[i] > 0)
                {
                    nDelete++;
                }
            }

            label npDeleteGlobal = nDelete;
            reduce(npDeleteGlobal, sumOp<label>());
            totalDeleteBeforeAdd_ += npDeleteGlobal;

            lammps_delete_particle(lmp_, deleteList, nDelete);

            delete [] deleteList;
        }

        addNewParticles();
        timeToAddParticle_ = addParticleTimeStep_;
    }
    else
    {
        timeToAddParticle_ -= runTime_.deltaT().value();
        Info<< "Time to add particle: " << timeToAddParticle_ << endl;
    }

    if (deleteParticleOption_ > 0)
    {
        int size = deleteParticleList_.size();
        int* deleteList = new int [size];
        int nDelete = 0;
        for (int i = 0; i < size; i++)
        {
            // TODO: this may have problem for parallel computing
            deleteList[i] = deleteParticleList_[i];
            if (deleteList[i] > 0)
            {
                nDelete++;
            }
        }

        label npDeleteGlobal = nDelete;
        reduce(npDeleteGlobal, sumOp<label>());
        totalDelete_ += npDeleteGlobal;

        lammps_delete_particle(lmp_, deleteList, nDelete);

        delete [] deleteList;
    }
}

//- find the CFD cell ID in which the particles are added
void softParticleCloud::findAddParticleCells()
{

    vector testPoint = vector::zero;
    Info << "test point: " << pointInRegion(testPoint, addParticleBox_) << endl;

    label nP = 0;
    forAll(mesh_.C(), cellI)
    {
        vector meshC = mesh_.C()[cellI];

        if (pointInRegion(meshC, addParticleBox_))
        {
            nP++;
        }
    }

    label nprocs = Pstream::nProcs();
    label myrank = Pstream::myProcNo();

    addParticleLocalList_.setSize(nprocs, 0);

    for(int i = 0; i < nprocs; i++)
    {
        addParticleLocalList_[i] = 0;
    }

    addParticleLocalList_[myrank] = nP;

    for(int i = 0; i < nprocs; i++)
    {
        reduce(addParticleLocalList_[i], maxOp<label>());
    }

    addParticleCellID_.setSize(nP);
    int i = 0;
    forAll(mesh_.C(), cellI)
    {
        vector meshC = mesh_.C()[cellI];

        if (pointInRegion(meshC, addParticleBox_))
        {
            addParticleCellID_[i] = cellI;
            i++;
        }
    }
}

bool softParticleCloud::pointInRegion(vector& point, tensor& box)
{
    scalar x1 = box.component(0);
    scalar x2 = box.component(1);
    scalar y1 = box.component(2);
    scalar y2 = box.component(3);
    scalar z1 = box.component(4);
    scalar z2 = box.component(5);
    scalar r1 = box.component(6);
    scalar r2 = box.component(7);

    if (addParticleOption_ == 1)
    {
        if
        (
            (point.x() - x1)*(point.x() - x2) < SMALL &&
            (point.y() - y1)*(point.y() - y2) < SMALL &&
            (point.z() - z1)*(point.z() - z2) < SMALL
        )
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else if (addParticleOption_ == 2)
    {
        vector p2p1 = vector(x2 - x1, y2 - y1, z2 - z1);
        scalar h = mag(p2p1);
        vector pxp1 = vector(point.x() - x1, point.y() - y1, point.z() - z1);
        scalar dot = (p2p1 & pxp1);

        vector p2p1E = p2p1;
        scalar hE = mag(p2p1E);
        vector pxp1E = pxp1 - addParticleBoxEccentricity_;
        scalar dotE = (p2p1E & pxp1E);

        if (dot < 0.0 || dot > pow(h,2))
        {
            return 0;
        }
        else
        {
            scalar dsq = (pxp1 & pxp1) - dot*dot/pow(h,2);
            scalar dsqE = (pxp1E & pxp1E) - dot*dot/pow(h,2);

            if (dsqE > r1*r1 && dsq < r2*r2)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

    }
    else
    {
        return 0;
    }
}


bool softParticleCloud::pointInBox(vector& point, tensor& box)
{
    scalar x1 = box.component(0);
    scalar x2 = box.component(1);
    scalar y1 = box.component(2);
    scalar y2 = box.component(3);
    scalar z1 = box.component(4);
    scalar z2 = box.component(5);

    if
    (
        (point.x() - x1)*(point.x() - x2) < SMALL &&
        (point.y() - y1)*(point.y() - y2) < SMALL &&
        (point.z() - z1)*(point.z() - z2) < SMALL
    )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void softParticleCloud::writeFields() const
{
    softParticle::writeFields(*this);
}

} // namespace Foam


// ************************************************************************* //
