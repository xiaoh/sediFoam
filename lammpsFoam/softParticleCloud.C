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
    nGlobal_ = lammps_get_natoms(lmp_);
    Info<< "FOAM reported # of particles according to Lammps: "
        << nGlobal_ << endl;

    // Setup temporary space for holding x & v for sending to Lammps.
    // TODO: find no. of particle on my processor
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
        label tags = scalar(tag[i]);
        label types = scalar(type[i]);

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
                    types
                );

            if (debug)
            {
                Pout<< "position is:" << pos << endl;
                Pout<< "cell is:" << cellI << endl;
                Pout<< "foam tag is:" << tags << endl;
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
void softParticleCloud::setPositionVelo(vector XLocal[], vector VLocal[])
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


// #define DEBUG_EVOLVE

// Call Lammps and evolve forward some steps.  Given particle
// positions, velocities, and drag.
// This is a wrapper for the Lammps functions which actually does
// the job.
void  softParticleCloud::lammpsEvolveForward
(
    vector* XLocal,
    vector* VLocal,
    vectorList FLocal,
    int nstep
)
{
    vectorList FGlobal(nGlobal_,vector::zero);

    int i = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();

        if (debug)
        {
            Pout<< "foam before evolve id is: " << p.ptag() << endl;
            Pout<< "foam before evolve type is: " << p.ptype() << endl;

            Pout<< "foam before evolve position is: " << p.position() << endl;
            Pout<< "foam before evolve velocity is: " << p.U() << endl;
        }

        // Here, the particle id in lammps starts from 1.
        // So, (p.ptag() - 1) would be the index in the FGlobal list.
        FGlobal[p.ptag() - 1] = vector::zero;
        FGlobal[p.ptag() - 1] = FLocal[i];

        if (debug)
        {
            Pout<< "FLocal[" << i << "] = "
                << FLocal[i] << endl;
            Pout<< "FGlobal[" << p.ptag() - 1 << "] = "
                << FGlobal[p.ptag() - 1] << endl;
            Pout<< "FGlobal[" << p.ptag() - 1 << "] position: "
                << p.position() << endl;
        }
    }

    Pstream::listCombineGather(FGlobal, plusEqOp<vector>());
    Pstream::listCombineScatter(FGlobal);

    flattenVectors(FGlobal, fArray_);

    if (debug)
    {
        for (label i = 0; i < nGlobal_; i++)
        {
                Pout<< "fArray_[" << 3*i + 0 << "] is: "
                    << fArray_[3*i + 0] << endl;
                Pout<< "fArray_[" << 3*i + 1 << "] is: "
                    << fArray_[3*i + 1] << endl;
                Pout<< "fArray_[" << 3*i + 2 << "] is: "
                    << fArray_[3*i + 2] << endl;
        }
    }

    lammps_put_drag(lmp_, fArray_); // Give current drag to Lammps

    // Ask lammps to move certain steps forward
    lammps_step(lmp_, nstep);

    // Harvest the lammps-updated x and v
    lammps_get_coord_velo(lmp_, xArray_, vArray_);

    double* xArrayLocal = new double [3*size()];
    double* vArrayLocal = new double [3*size()];

    i = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();

        xArrayLocal[3*i] = xArray_[3*p.ptag() - 3];
        xArrayLocal[3*i + 1] = xArray_[3*p.ptag() - 2];
        xArrayLocal[3*i + 2] = xArray_[3*p.ptag() - 1];

        vArrayLocal[3*i] = vArray_[3*p.ptag() - 3];
        vArrayLocal[3*i + 1] = vArray_[3*p.ptag() - 2];
        vArrayLocal[3*i + 2] = vArray_[3*p.ptag() - 1];

        if (debug)
        {
            Pout<< "position of particle# " << i << " is: "
                << xArrayLocal[3*i] << " " << xArrayLocal[3*i + 1] << " "
                << xArrayLocal[3*i + 2] << endl;
        }
    }

    // This a temperary work-around. Assemble array to vectors.
    assembleVectors(XLocal, xArrayLocal);
    assembleVectors(VLocal, vArrayLocal);

    delete [] xArrayLocal;
    delete [] vArrayLocal;

} // Job done; Proceed to next fluid calculation step.


void softParticleCloud::writeFields() const
{
    softParticle::writeFields(*this);
}

} // namespace Foam


// ************************************************************************* //
