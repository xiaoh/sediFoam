/*---------------------------------------------------------------------------* \
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

#include "enhancedCloud.H"
// #define DEBUG_FORCE


namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Setup particle diameter list
//  called in constructor
void enhancedCloud::setupParticleDia()
{
    pDia_.setSize(particleCount_);
    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();
        pDia_[particleI] = p.d();
    }
}


//- Update particle alpha list (per fluid step)
void  enhancedCloud::updateParticleAlpha()
{
    pAlpha_.setSize(particleCount_);
    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();
        label cellI = p.cell();
        if (cellI<0)
        {
            Info << "cell not found!" << endl;
            continue;
        }
        pAlpha_[particleI] = gamma_[cellI];
    }
}


//- Update relative Ur list (per fluid step or substep)
//  Update Ur_ and magUr_
//  Called by calcTcField before calculating Omega & Asrc
//  Also called after each substep.
void enhancedCloud::updateParticleUr()
{
    Uri_.setSize(particleCount_);
    magUri_.setSize(particleCount_);

    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();
        if (p.cell() < 0)
        {
            Uri_[particleI] = vector::zero;
            magUri_[particleI] = scalar(0);
            ++pIter;
            continue;
        }

        // velocity of the fluid
        Uri_[particleI] = UfSmoothed_[p.cell()] - p.U();
        magUri_[particleI] = mag(Uri_[particleI]);
    }
}


void  enhancedCloud::updateDragOnParticles()
{
    vectorField gradp = fvc::grad(pf_)().internalField();
    vectorField curlU = fvc::curl(Uf_)().internalField();

    setupParticleDia();
    updateParticleAlpha();

    pDrag_.setSize(particleCount_);
    pDuDt_.setSize(particleCount_);
    Jd_.setSize(particleCount_);

    // clear pDrag_
    pDrag_ = vector::zero;
    pDuDt_ = vector::zero;
    Jd_ *= 0.0;

    Jd_ = drag_->Jd(magUri_);

    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();

        if (p.cell() < 0) { continue;}

        if (debug)
        {
            Pout<< "---Jd ";
            Pout<< Jd_[particleI];
            Pout<< "  particle#: " << particleI
                << "  Cell I: " << p.cell()
                << "  alpha " << pAlpha_[particleI]
                << "  Uri_: "  << Uri_[particleI] << endl;
        }

        // local pressure gradient is used though (no weighting).
        pDrag_[particleI] = vector::zero;
        pDuDt_[particleI] = DDtUf_[p.cell()];

        if (particleDragFlag_)
        {
            pDrag_[particleI] +=
                Jd_[particleI]*(1.0 - pAlpha_[particleI])
               *p.Vol()*Uri_[particleI];        // Drag
        }
        if (particlePressureGradFlag_)
        {
            vector gradpVec = gradp[p.cell()];

            pDrag_[particleI] += -gradpVec*p.Vol();             // Dynamic pressure
        }
        if (particleBuoyancyFlag_)
        {
            pDrag_[particleI] += 
               -gravity_*rhob_*p.Vol();      // Buoyancy
        }
        // Added mass force
        if (particleAddedMassFlag_)
        {
            vector dupdt = (p.U()-p.UOld())/runTime().deltaT().value();
            vector acceleration = DDtUf_[p.cell()]-dupdt;
            if (mag(DDtUf_[p.cell()] - dupdt) > 10)
            {
                // Avoid too large added mass
                // dupdt = dupdt/mag(dupdt)*100;
                acceleration =  acceleration/(mag(acceleration)+ROOTVSMALL)*10;
            }

            // Add mass is only effective when the particle is accelerating
            pDrag_[particleI] += 0.5*rhob_*p.Vol()*acceleration;
        }
        if (particleLiftForceFlag_)
        {
            scalar liftCoeff = 1.6;
            pDrag_[particleI] +=
                liftCoeff*rhob_*sqrt(nub_)*sqr(p.d())
             *((Uri_[particleI])^curlU[p.cell()])
               /sqrt(mag(curlU[p.cell()]) + ROOTVSMALL);
        }
        if (particleHistoryForceFlag_)
        {
        //    unfinished
        }
        if (lubricationFlag_)
        {
            scalar distMin = 0.0001*p.d();
            scalar distMax = 0.1*p.d();
            scalar distWall = p.position().y() - 0.5*p.d();
            scalar pVel = p.U().y();
            if (distWall < distMax && distWall > distMin)
            {
                vector normalVec = vector(0,1,0);
                pDrag_[particleI] +=
                    6*3.1416*nub_*rhob_
                   *(-pVel)/distWall*(p.d()*p.d())/4.0*normalVec;
            }
        }
        if (mag(inletForceRatio_) > 0)
        {
            if (pointInRegion(p.position(), inletBox_))
            {
                pDrag_[particleI] =
                    p.m()*(inletForceRatio_-p.U())
                   /runTime().deltaT().value();
            }
        }

#ifdef DEBUG_FORCE
        Info<< "drag is: "
            << Jd_[particleI]*(1.0 - pAlpha_[particleI])
              *p.Vol()*Uri_[particleI] << endl;
        Info<< "pressure gradient force is: "
            << - gradp[p.cell()]*p.Vol() << endl;

        vector dupdt = (p.U()-p.UOld())/runTime().deltaT().value();
        Info<< "volume of particle is: " << p.Vol() << endl;
        Info<< "current velocity is: " << p.U() << endl;
        Info<< "previous velocity is: " << p.UOld() << endl;
        Info<< "current acceleration is: " << dupdt << endl;
        if (mag(dupdt) > 100)
        {
            dupdt = dupdt/mag(dupdt)*100;
            Info<< "acceleration too large" << endl;
            Info<< "the value is adjusted to: " << dupdt << endl;
        }
        Info<< "added mass is: "
            << - 0.5*rhob_*p.Vol()*dupdt<< endl;
        Info<< "total force is: " << pDrag_[particleI] << endl;


        if (particleLiftForceFlag_)
        {
            Info<< "Uri is: " << Uri_[particleI] << endl;
            Info<< "curlU is: " << curlU[p.cell()] << endl;
            Info<< "Uri x curlU is: "
                << (Uri_[particleI]^curlU[p.cell()]) << endl;
            Info<< "d^2 is: " << sqr(p.d()) << endl;

            Info<< "lift force is: "
                <<  1.6*rhob_*sqrt(nub_)*sqr(p.d())
                   *((Uri_[particleI])^curlU[p.cell()])
                   /sqrt(mag(curlU[p.cell()])) << endl;
        }
        if (lubricationFlag_)
        {
            scalar distMin = 0.0001*p.d();
            scalar distMax = 0.1*p.d();
            scalar distWall = p.position().y() - 0.5*p.d();
            scalar pVel = p.U().y();
            if (distWall < distMax && distWall > distMin)
            {
                vector normalVec = vector(0,1,0);
                Info<< "lubrication is: "
                    << 6*3.1416*nub_*rhob_
                      *(-pVel)/distWall*(p.d()*p.d())/4.0*normalVec
                    << endl;
            }
        }
#endif
    }
}


//- Use averaging method to get Omega and Asrc field
void enhancedCloud::calcTcFields()
{
    Omega_.internalField() *= 0.0;
    Asrc_.internalField() *= 0.0;
    Asrc2_.internalField() *= 0.0;

    // prepare alpha & Ur list for drag computation
    updateParticleAlpha();
    updateParticleUr();
    setupParticleDia();

    if (debug)
    {
        Info<< "particleAlpha Field updated: "
            << pAlpha_ << endl;
        Info<< "particleUr Field updated: "
            << Uri_ << " mag: "
            << magUri_ << endl;
    }

    Jd_ = drag_->Jd(magUri_);

    bool semiImplicit = 0;
    if (semiImplicit)
    {
        // Scan particle list to average and get Omega & Asrc field
        label particleI = 0;
        for
        (
            softParticleCloud::iterator pIter = softParticleCloud::begin();
            pIter != softParticleCloud::end();
            ++pIter, ++particleI
        )
        {
            softParticle& p = pIter();

            // Update Omega field (fluid density omitted)
            label cellI = p.cell();
            if(cellI < 0) continue;
            scalar omg = p.Vol()*Jd_[particleI]/(mesh_.V()[cellI]);
            Omega_.internalField()[cellI] += omg;
            Asrc_.internalField()[cellI] += omg * p.U();
            ++pIter;
        }
    }
    else
    {
        // scan particle list to average and get Omega & Asrc field
        label particleI = 0;
        for
        (
            softParticleCloud::iterator pIter = softParticleCloud::begin();
            pIter != softParticleCloud::end();
            ++pIter, ++particleI
        )
        {
            softParticle& p = pIter();

            // update Omega field (fluid density omitted)
            label cellI = p.cell();

            if (cellI < 0) continue;

            scalar omg = p.Vol()*Jd_[particleI]/(mesh_.V()[cellI]);

            // accumulate drag from particles to host cells
            // to be smoothed later!
            Omega_.internalField()[cellI] += omg;
            Asrc_.internalField()[cellI] += omg*(p.U() - UfSmoothed_[cellI]);
                                          // + gravity_*rhob_*p.Vol()/(1-gamma_[cellI])/(mesh_.V()[cellI]);
            Asrc2_.internalField()[cellI] += omg*(p.U() - UfSmoothed_[cellI]);
                                          // + gravity_*rhob_*p.Vol()/(1-gamma_[cellI])/(mesh_.V()[cellI]);
            // Asrc_.internalField()[cellI] += omg*(p.U() - Uf_[cellI]);
        }

        Omega_.internalField() *= 0;

        // F1 and F2 are calculated to show that
        // the momentum is conservative
        vector Ftotal1(vector::zero);

        forAll(Asrc_.internalField(), ceI)
        {
            Ftotal1 +=
                Asrc_.internalField()[ceI]*mesh_.V()[ceI]
              *(1 - gamma_.internalField()[ceI]);
        }

        // TODO: Derive the conservative way of diffusion;
        // implement, and check numerically.
        // Smoothing operations
        Asrc_.internalField() =
                Asrc_.internalField()*(1 - gamma_.internalField());

        if (dragSmoothFlag_)
        {
            smoothField(Asrc_);
        }

        Asrc_.internalField() /=
            (1 - gamma_.internalField());

        // TODO: assert particle total force equal Eulerian field total force

        // F1 and F2 are calculated to show that
        // the momentum is conservative
        vector Ftotal2(vector::zero);

        forAll(Asrc_.internalField(), ceI)
        {
            Ftotal2 +=
                Asrc_.internalField()[ceI]*mesh_.V()[ceI]
              *(1 - gamma_.internalField()[ceI]);
        }

        reduce(Ftotal1, sumOp<vector>());
        reduce(Ftotal2, sumOp<vector>());

        Info<< "total F before: " << Ftotal1 << endl;
        Info<< "total F after: " << Ftotal2 << endl;

        Asrc_.correctBoundaryConditions();
        Omega_.correctBoundaryConditions();
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from component
enhancedCloud::enhancedCloud
(
    const volVectorField& U,
    const volScalarField& p,
    volVectorField& Ue,
    const volVectorField& Uf,
    const volVectorField& DDtUf,
    dimensionedScalar nu,
    volScalarField& alpha,
    IOdictionary& cloudDict,
    IOdictionary& transDict,
    scalar diffusionBandWidth,
    label diffusionSteps
)
:
    softParticleCloud(U, p, Ue, nu, alpha, cloudDict),
    mesh_(U.mesh()),
    Uf_(Uf),
    DDtUf_(DDtUf),
    UfSmoothed_(Uf),
    Omega_
    (
        IOobject
        (
            "Omega",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0),
            scalar(0.0)
        )
    ),
    Asrc_
    (
        IOobject
        (
            "A_Source",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
           "zero",
           dimensionSet(1, -2, -2, 0, 0),
           vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    ),
    Asrc2_
    (
        IOobject
        (
            "A_Source_unsmoothed",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
           "zero",
           dimensionSet(1, -2, -2, 0, 0),
           vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    ),
    diffusionRunTime_
    (
        "controlDiffDict",
        runTime().rootPath(),
        runTime().caseName()
    ),
    diffusionMesh_
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            diffusionRunTime_.timeName(),
            diffusionRunTime_,
            Foam::IOobject::MUST_READ
        )
    ),
    simple_(diffusionMesh_),
    diffusionTimeCount_(2, 0.0),
    particleMoveTime_(0.0)
{
    drag_ = Foam::dragModel::New(cloudDict, transDict, pAlpha_, pDia_);
    dimensionedScalar rhob(transDict.lookup("rhob"));
    rhob_ = rhob.value();

    dimensionedScalar nub(transDict.lookup("nub"));
    nub_ = nub.value();

    particleCount_ = size();

    forAllIter(softParticleCloud, *this, iter)
    {
        softParticle& p = iter();

        label cellI = p.cell();

        // alpha field
        gamma_.internalField()[cellI] +=
            p.Vol();
    }

    // determine the time and time step in diffusion procedure
    scalar diffusionTime = pow(diffusionBandWidth, 2)/4;
    scalar diffusionDeltaT = diffusionTime/(diffusionSteps + ROOTVSMALL);

    diffusionRunTime_.setEndTime(diffusionTime);
    diffusionRunTime_.setDeltaT(diffusionDeltaT);
    Info << "diffusion time is: " << diffusionTime << endl;
    Info << "diffusion time step is: " << diffusionDeltaT << endl;

    // determine the fields to diffuse
    UfSmoothFlag_ = cloudProperties_.lookupOrDefault("UfSmooth", true);
    UpSmoothFlag_ = cloudProperties_.lookupOrDefault("UpSmooth", true);
    dragSmoothFlag_ = cloudProperties_.lookupOrDefault("dragSmooth", true);
    alphaSmoothFlag_ = cloudProperties_.lookupOrDefault("alphaSmooth", true);

    smoothDirection_ =
        cloudProperties_.lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        );

    // determine the forces to add
    particleDragFlag_ = cloudProperties_.lookupOrDefault("particleDrag", true);
    particlePressureGradFlag_ =
        cloudProperties_.lookupOrDefault("particlePressureGrad", true);
    particleBuoyancyFlag_ =
        cloudProperties_.lookupOrDefault("particleBuoyancy", false);
    particleAddedMassFlag_ =
        cloudProperties_.lookupOrDefault("particleAddedMass", false);
    particleLiftForceFlag_ =
        cloudProperties_.lookupOrDefault("particleLift", false);
    particleHistoryForceFlag_ =
        cloudProperties_.lookupOrDefault("particleHistoryForce", false);
    lubricationFlag_ =
        cloudProperties_.lookupOrDefault("lubricationForce", false);

    if (addParticleOption_ > 0)
    {
        inletForceRatio_ =
            cloudProperties_.lookupOrDefault("inletForce", vector::zero);
    }
    else
    {
        inletForceRatio_ = vector::zero;
    }


    Info<< "forces are: "
        << particleDragFlag_
        << particlePressureGradFlag_
        << particleBuoyancyFlag_
        << particleAddedMassFlag_
        << particleLiftForceFlag_
        << particleHistoryForceFlag_
        << lubricationFlag_
        << inletForceRatio_
        << endl;

    // initial quantities for dragModel
    pDia_.setSize(particleCount_);
    pAlpha_.setSize(particleCount_);
    Uri_.setSize(particleCount_);
    magUri_.setSize(particleCount_);

    // initialise drag force on each particle
    pDrag_ = vectorList(particleCount_);

    // initialise DuDt of flow velocity on each particle
    pDuDt_ = vectorList(particleCount_);

    // initialize alpha and Ue field
    particleToEulerianField();

    setupParticleDia();

    // smooth fluid velocity to initialise the lift&drag coefficients
    // gamma is smoothed field since we smooth it in particleToEulerianField();
    UfSmoothed_.internalField() = Uf_.internalField();

    if (UfSmoothFlag_)
    {
        UfSmoothed_.internalField() *=
            (1 - gamma_.internalField());

        Info<< "smooth Uf..." << endl;
        smoothField(UfSmoothed_);

        UfSmoothed_.internalField() /=
            (1 - gamma_.internalField());

        UfSmoothed_.correctBoundaryConditions();
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

enhancedCloud::~enhancedCloud()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void enhancedCloud::evolve()
{
    softParticle::trackingData td0(*this);

    label Ns = subCycles_;

    UfSmoothed_.internalField() = Uf_.internalField();

    // smooth Uf when necessary
    if (UfSmoothFlag_)
    {
        UfSmoothed_.internalField() *=
            (1 - gamma_.internalField());

        Info<< "smooth Uf..." << endl;
        smoothField(UfSmoothed_);

        UfSmoothed_.internalField() /=
            (1 - gamma_.internalField());

        UfSmoothed_.correctBoundaryConditions();
    }


    // evolve Ns steps forward each time when Lammps is called.
    for (label k = 0; k < Ns; k++)
    {

        // adding particles when the conditions are satisfied
        if (addParticleOption_ > 0 && timeToAddParticle_ <= 0)
        {
            if (deleteBeforeAddFlag_)
            {
                deleteParticleBeforeAdd();
            }

            addParticleOpenFOAM();
        }
        // deleting particles when the conditions are satisfied
        if (deleteParticleOption_ > 0)
        {
            deleteParticleOpenFOAM();
        }

        label nLocal = size(); // Get number of local particles

        if (particleCount_ != size())
        {
            particleCount_ = size();
        }

        // the squence of the data on XLocal and VLocal is the same as
        // sequence of local particle index.
        vector* XLocal = new vector[nLocal];
        vector* VLocal = new vector[nLocal];
        int* lmpCpuIdLocal = new int[nLocal];

        int nstep = subSteps_;

        // here, Ur (relative velocity is called)
        updateParticleUr();

        updateDragOnParticles();

        // XLocal & VLocal are work spaces for "lammpsEvolveForward"
        // newly obtianed values are put there
        lammpsEvolveForward
        (
            XLocal,
            VLocal,
            lmpCpuIdLocal,
            pDrag_,
            pDuDt_,
            nstep
        );

        // update position/velocity of all particles in this cloud.
        // (Harvest XLocal & VLocal)  Lammps --> Cloud
        setPositionVeloCpuId(XLocal, VLocal, lmpCpuIdLocal);

        diffusionRunTime_.cpuTimeIncrement();
        // move particle to the new position
        Cloud<softParticle>::move(td0, mesh_.time().deltaTValue());

        particleMoveTime_ += diffusionRunTime_.cpuTimeIncrement();

        if (particleCount_ != size())
        {
            // Pout<< "Warning: enhancedCloud::evolve: "
            //     << "particle number modified! "
            //     << "Particle number before: " << particleCount_
            //     << "Particle number now: " << size()
            //     << endl;

            particleCount_ = size();
        }

        // make sure all particles are in cell.
        // assertParticleInCell();

        // update Uri (relative particle velocities) here.
        updateParticleUr();

        // change Eulerian (mesh-based) alpha field
        if (k == 0)
        {
            particleToEulerianField();
        }

        delete [] XLocal;
        delete [] VLocal;
        delete [] lmpCpuIdLocal;
    }

    // Pout<< "After this cycle, "
    //     << size() << " local particles has been moved. " << endl;
    Info<< "Adding/deleting statistics. adding: " << totalAdd_ << " deleting: " << totalDelete_
        << " deleting before add: " << totalDeleteBeforeAdd_ << endl;
}


void enhancedCloud::smoothField(volScalarField& sFieldIn)
{
    scalar t0 = runTime().elapsedCpuTime();
    volScalarField diffWorkField
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            dimless,
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
    );

    diffWorkField.internalField() = sFieldIn.internalField();

    dimensionedTensor DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);

    scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

    Info<< "smoothing " << sFieldIn.name() << endl;

    diffusionTimeCount_[0] += runTime().elapsedCpuTime() - t0;
    t0 = runTime().elapsedCpuTime();

    while (diffusionRunTime_.loop())
    {
        if (diffusionRunTime_.timeIndex() == 1)
        {
            while (simple_.correctNonOrthogonal())
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }
        else
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }

    diffusionTimeCount_[1] += runTime().elapsedCpuTime() - t0;
    t0 = runTime().elapsedCpuTime();

    diffusionRunTime_.setTime(startTime,startIndex);

    sFieldIn.internalField() = diffWorkField.internalField();
    diffusionTimeCount_[0] += runTime().elapsedCpuTime() - t0;
}


void enhancedCloud::smoothField(volVectorField& sFieldIn)
{
    scalar t0 = runTime().elapsedCpuTime();
    volVectorField diffWorkField
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            dimless,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );

    diffWorkField.internalField() = sFieldIn.internalField();

    dimensionedTensor DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_);

    scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

    Info<< "smoothing " << sFieldIn.name() << endl;

    diffusionTimeCount_[0] += runTime().elapsedCpuTime() - t0;
    t0 = runTime().elapsedCpuTime();

    while (diffusionRunTime_.loop())
    {
        if (diffusionRunTime_.timeIndex() == 1)
        {
            while (simple_.correctNonOrthogonal())
            {
                solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
            }
        }
        else
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }

    diffusionTimeCount_[1] += runTime().elapsedCpuTime() - t0;
    t0 = runTime().elapsedCpuTime();

    diffusionRunTime_.setTime(startTime,startIndex);

    sFieldIn.internalField() = diffWorkField.internalField();
    diffusionTimeCount_[0] += runTime().elapsedCpuTime() - t0;
}


//- Refresh Ue and Gamma using Gaussian averaging
void enhancedCloud::particleToEulerianField()
{
    Info<< "particle to eulerian field ..." << endl;
    gamma_.internalField() *= 0.0;
    Ue_.internalField() *= 0.0;

    // obtain alpha and Ua field
    forAllIter(softParticleCloud, *this, iter)
    {
        softParticle& p = iter();

        label cellI = p.cell();

        // alpha field
        gamma_.internalField()[cellI] += p.Vol();

        Ue_.internalField()[cellI] += p.Vol()*p.U();
    }

    gamma_.internalField() /= mesh_.V();

    // Utotal1 and Utotal2 are calculated to show that
    // the momentum is conservative
    vector Utotal1(vector::zero);

    forAll(Ue_.internalField(), ceI)
    {
        Utotal1 += Ue_.internalField()[ceI];
    }

    Ue_.internalField() /= mesh_.V();

    // smooth alpha and Ua field
    if (alphaSmoothFlag_)
    {
        Info<< "smoothing alpha flag on..." << endl;
        smoothField(gamma_);
    }

    if (UpSmoothFlag_)
    {
        smoothField(Ue_);
    }

    forAll(Ue_.internalField(), ceI)
    {
        if (gamma_.internalField()[ceI] > ROOTVSMALL)
        {
            Ue_.internalField()[ceI] /=
                (gamma_.internalField()[ceI]);
        }
    }

    vector Utotal2(vector::zero);

    forAll(Ue_.internalField(), ceI)
    {
        Utotal2 +=
            Ue_.internalField()[ceI]*mesh_.V()[ceI]*gamma_.internalField()[ceI];
    }

    reduce(Utotal1, sumOp<vector>());
    reduce(Utotal2, sumOp<vector>());

    Info<< "total U solid before: " << Utotal1 << endl;
    Info<< "total U solid after: " << Utotal2 << endl;

    Ue_.correctBoundaryConditions();
    gamma_.correctBoundaryConditions();
}


//- Check if all particles in fluid cells
//  It is in this class because it used to depend on neighbour info.
//  Not any more though.
void enhancedCloud::assertParticleInCell()
{
    softParticleCloud::iterator pIter = begin();

    forAllIter(softParticleCloud, *this, iter)
    {
        softParticle& p = pIter();
        point& pos = p.position();

        bool found = false;

        label lCellI = mesh_.findCell(pos);

        if (lCellI >= 0) found = true;

        if (!found)
        {
            WarningIn("enhancedCloud::assertParticleInCell()")
                << "Particle outside fluid domain!" << endl
                << " Debug Info: "
                << " Particle Position: " << pos
                << " Associated cell ID: " << p.cell()
                << abort(FatalError);
        }
    }

}

//- Add OpenFOAM particles
void enhancedCloud::addParticleOpenFOAM()
{
    maxTag_ = 0;
    int i = 0;
    for
    (
        softParticleCloud::iterator pIter = begin();
        pIter != end();
        ++pIter, ++i
    )
    {
        softParticle& p = pIter();
        maxTag_ = max(p.ptag(),maxTag_);
    }

    reduce(maxTag_, maxOp<label>());

    label myrank = Pstream::myProcNo();

    forAll(addParticleLocalList_, i)
    {
        if (myrank > i)
        {
            maxTag_ += addParticleLocalList_[i];
        }
    }

    forAll(addParticleCellID_, i)
    {

        label cellI = addParticleCellID_[i];
        vector pos = mesh_.C()[cellI];

        vector velo = addParticleVel_;

        scalar ds = addParticleInfo_[0];
        scalar rhos = addParticleInfo_[1];
        label types = int(addParticleInfo_[2]);

        // TODO: this may have problems for parallel computing
        label lmpCpuIds = 0;
        forAll(lmpLocalBoxList_, boxI)
        {
            if (pointInBox(pos, lmpLocalBoxList_[boxI]))
            {
                lmpCpuIds = boxI;
            }
        }

        label tags = maxTag_ + 1 + i;

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

        softParticleCloud::addParticle(ptr);
    }
}


//- Delete OpenFOAM particles
void enhancedCloud::deleteParticleBeforeAdd()
{
    if (deleteParticleOption_ > 0)
    {

        label nprocs = Pstream::nProcs();

        int i = 0;
        int nDelete = 0;
        for
        (
            softParticleCloud::iterator pIter = begin();
            pIter != end();
            ++pIter, ++i
        )
        {
            softParticle& p = pIter();
            vector pPosition = p.position();

            if (pointInRegion(pPosition, clearInitialBox_) && (p.ptype() == 1))
            {
                nDelete++;
            }
        }

        deleteBeforeAddList_.setSize(nDelete);
        List<labelList> deleteBeforeAddListList(nprocs);
        labelList deletedParticleNo(nprocs, 0);
        labelList assembleLmpCpuIdList(nDelete, 0);

        for (int i = 0; i < nDelete; i++)
        {
            // TODO: this may have problem for parallel computing
            deleteBeforeAddList_[i] = 0;
        }

        i = 0;
        for
        (
            softParticleCloud::iterator pIter = begin();
            pIter != end();
            ++pIter
        )
        {
            softParticle& p = pIter();
            vector pPosition = p.position();

            if (pointInRegion(pPosition, clearInitialBox_) && (p.ptype() == 1))
            {
                Info<< "deleting particle.." << endl;
                // TODO: this may have problem for parallel computing
                deleteBeforeAddList_[i] = p.ptag();
                label boxI = p.pLmpCpuId();
                deleteParticle(p);
                deletedParticleNo[boxI] ++;
                assembleLmpCpuIdList[i] = boxI;
                i++;
            }
        }

        assembleList<labelList>
        (
            deleteBeforeAddList_,
            deleteBeforeAddListList,
            deletedParticleNo,
            assembleLmpCpuIdList
        );

        List<labelList> toLmpDeleteTagListList(nprocs);

        transposeAmongProcs<labelList>
        (
            deleteBeforeAddListList,
            toLmpDeleteTagListList
        );

        label toLmpListSize = 0;
        forAll(deleteBeforeAddListList, listI)
        {
            toLmpListSize += toLmpDeleteTagListList[listI].size();
        }
        labelList toLmpDeleteTagList(toLmpListSize, 0);

        // flattenList<labelList> (toLmpDeleteTagListList, toLmpDeleteTagList);

        i = 0;
        forAll(toLmpDeleteTagListList,listI)
        {
            forAll(toLmpDeleteTagListList[listI], memberI)
            {
                toLmpDeleteTagList[i] = toLmpDeleteTagListList[listI][memberI];
                i++;
            }
        }

        deleteBeforeAddList_.setSize(toLmpListSize);

        forAll(deleteBeforeAddList_, i)
        {
            deleteBeforeAddList_[i] = toLmpDeleteTagList[i];
        }
    }
}


//- Delete OpenFOAM particles
void enhancedCloud::deleteParticleOpenFOAM()
{
    if (deleteParticleOption_ > 0)
    {
        label nprocs = Pstream::nProcs();

        int i = 0;
        int nDelete = 0;
        for
        (
            softParticleCloud::iterator pIter = begin();
            pIter != end();
            ++pIter, ++i
        )
        {
            softParticle& p = pIter();
            vector pPosition = p.position();

            if (pointInRegion(pPosition, deleteParticleBox_) && (p.ptype() == 1))
            {
                nDelete++;
            }
        }

        deleteParticleList_.setSize(nDelete);
        List<labelList> deleteParticleListList(nprocs);
        labelList deletedParticleNo(nprocs, 0);
        labelList assembleLmpCpuIdList(nDelete, 0);

        for (int i = 0; i < nDelete; i++)
        {
            // TODO: this may have problem for parallel computing
            deleteParticleList_[i] = 0;
        }

        i = 0;
        for
        (
            softParticleCloud::iterator pIter = begin();
            pIter != end();
            ++pIter
        )
        {
            softParticle& p = pIter();
            vector pPosition = p.position();

            if (pointInRegion(pPosition, deleteParticleBox_) && (p.ptype() == 1))
            {
                Info<< "deleting particle.." << endl;
                // TODO: this may have problem for parallel computing
                deleteParticleList_[i] = p.ptag();
                label boxI = p.pLmpCpuId();
                deleteParticle(p);
                deletedParticleNo[boxI] ++;
                assembleLmpCpuIdList[i] = boxI;
                i++;
            }
        }

        assembleList<labelList>
        (
            deleteParticleList_,
            deleteParticleListList,
            deletedParticleNo,
            assembleLmpCpuIdList
        );

        List<labelList> toLmpDeleteTagListList(nprocs);

        transposeAmongProcs<labelList>
        (
            deleteParticleListList,
            toLmpDeleteTagListList
        );

        label toLmpListSize = 0;
        forAll(deleteParticleListList, listI)
        {
            toLmpListSize += toLmpDeleteTagListList[listI].size();
        }

        labelList toLmpDeleteTagList(toLmpListSize, 0);

        // flattenList<labelList> (toLmpDeleteTagListList, toLmpDeleteTagList);

        i = 0;
        forAll(toLmpDeleteTagListList,listI)
        {
            forAll(toLmpDeleteTagListList[listI],memberI)
            {
                toLmpDeleteTagList[i] = toLmpDeleteTagListList[listI][memberI];
                i++;
            }
        }

        deleteParticleList_.setSize(toLmpListSize);

        forAll(deleteParticleList_, i)
        {
            deleteParticleList_[i] = toLmpDeleteTagList[i];
        }
    }
}


//- Display sum of drag as Eulerian field
void enhancedCloud::dragInfo()
{
    volVectorField dragSumField_
    (
        IOobject
        (
            "Drag_Sum_up",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 1, -2, 0, 0),
            vector::zero
        )
    );

    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();

        if (p.cell() < 0) continue;
        dragSumField_.internalField()[p.cell()] += pDrag_[particleI];
    }

    // convert to force on unit cell volume
    dragSumField_.internalField() /= mesh_.V();

    Info<< "Sum up drag on particle in unit cell: "
        << dragSumField_.internalField()
        << endl;
}

void enhancedCloud::averageInfo()
{
    vector averageVel(vector::zero);
    vector totalVel(vector::zero);
    scalar localVolume(0);
    scalar totalVolume(0);

    label particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();

        localVolume += p.Vol();
        totalVel += p.U()*p.Vol();
    }

    reduce(totalVel, sumOp<vector>());
    reduce(localVolume, sumOp<scalar>());
    totalVolume = localVolume;
    averageVel = totalVel/(totalVolume+ROOTVSMALL);

    Info<< "total volume of particles is: " << totalVolume << endl;
    Info<< "total (velocity x volume) of particles is: " << totalVel << endl;
    Info<< "average velocity of all particles is: " << averageVel << endl;
}


} // End namespace Foam


// ************************************************************************* //
