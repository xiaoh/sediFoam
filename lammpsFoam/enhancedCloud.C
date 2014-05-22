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
            Pout << "cell not found!" << endl;
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

    setupParticleDia();
    updateParticleAlpha();

    pDrag_.setSize(particleCount_);
    Jd_.setSize(particleCount_);

    // clear pDrag_
    pDrag_ = vector::zero;
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
        pDrag_[particleI] =
            Jd_[particleI]*(1.0 - pAlpha_[particleI])
           *p.Vol()*Uri_[particleI]         // Drag
          - gradp[p.cell()]*p.Vol();        // Buoyancy
    }
}


//- Use averaging method to get Omega and Asrc field
void enhancedCloud::calcTcFields()
{
    Omega_.internalField() *= 0.0;
    Asrc_.internalField() *= 0.0;

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
            Asrc_.internalField()[cellI] += omg*(p.U() - UfSmoothed_[cellI]);
        }

        Omega_.internalField() *= 0;

        // F1 and F2 are calculated to show that
        // the momentum is conservative
        vector Ftotal1(vector::zero);

        forAll(Asrc_.internalField(), ceI)
        {
            Ftotal1 += 
                Asrc_.internalField()[ceI]*mesh_.V()[ceI]*(1 - gamma_.internalField()[ceI]);
        }

        // TODO: Derive the conservative way of diffusion; implement, and check numerically.
        // Smoothing operations
        Asrc_.internalField() =
                Asrc_.internalField()*(1 - gamma_.internalField());
        smoothField(Asrc_);

        Asrc_.internalField() /=
            (1 - gamma_.internalField());

        // TODO: assert particle total force equal Eulerian field total force

        // F1 and F2 are calculated to show that
        // the momentum is conservative
        vector Ftotal2(vector::zero);

        forAll(Asrc_.internalField(), ceI)
        {
            Ftotal2 += 
                Asrc_.internalField()[ceI]*mesh_.V()[ceI]*(1 - gamma_.internalField()[ceI]);
        }

        reduce(Ftotal1, sumOp<vector>());
        reduce(Ftotal2, sumOp<vector>());

        Info<< "total F before: " << Ftotal1 << endl;
        Info<< "total F after: " << Ftotal2 << endl;
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from component
enhancedCloud::enhancedCloud
(
    const volPointInterpolation& vpi,
    const volVectorField& U,
    const volScalarField& p,
    volVectorField& Ue,
    const volVectorField& Uf,
    dimensionedScalar nu,
    volScalarField& alpha,
    IOdictionary& cloudDict,
    IOdictionary& transDict,
    scalar bwDxRatio,
    scalar diffusionBandWidth,
    scalar diffusionSteps
)
:
    softParticleCloud(vpi, U, p, Ue, nu, alpha, cloudDict),
    mesh_(U.mesh()),
    Uf_(Uf),
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
    diffusionTimeCount_(0.0),
    particleMoveTime_(0.0)
{
    drag_ = Foam::dragModel::New(cloudDict, transDict, pAlpha_, pDia_);

    particleCount_ = size();

    forAllIter(softParticleCloud, *this, iter)
    {
        softParticle& p = iter();

        label cellI = p.cell();

        // alpha field
        gamma_.internalField()[cellI] +=
            p.Vol();
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        vector(1,1,1) & mesh_.Sf()
    );

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().internalField()
    );

    scalar courantTimeStep =
        1/(0.5*pow(gMax(sumPhi/mesh_.V().field()),2));

    Info<< "Best time step is: " << courantTimeStep << endl;

    scalar diffusionTime = pow(diffusionBandWidth,2)/4;

    scalar diffusionDeltaT = courantTimeStep;
    scalar nDiffusionTimeStep = ceil(diffusionTime/diffusionDeltaT);
    Info<< "explicit DiffusionTimeStep is: " << nDiffusionTimeStep << endl;
    Info<< "implicit DiffusionTimeStep is: " << diffusionSteps << endl;

    diffusionDeltaT = diffusionTime/(diffusionSteps + SMALL);
    diffusionRunTime_.setEndTime(diffusionTime);
    diffusionRunTime_.setDeltaT(diffusionDeltaT);
    Info<< "diffusion time is: " << diffusionTime << endl;
    Info<< "diffusion time step is: " << diffusionDeltaT << endl;

    // initial quantities for dragModel
    pDia_.setSize(particleCount_);
    pAlpha_.setSize(particleCount_);
    Uri_.setSize(particleCount_);
    magUri_.setSize(particleCount_);

    // initialise drag force on each particle
    pDrag_ = vectorList(particleCount_);


    // initialize alpha and Ue field
    particleToEulerianField();

    setupParticleDia();

    // smooth fluid velocity to initialise the lift&drag coefficients
    // gamma is smoothed field since we smooth it in particleToEulerianField();
    UfSmoothed_.internalField() =
        Uf_.internalField()*(1 - gamma_.internalField());
    smoothField(UfSmoothed_);

    UfSmoothed_.internalField() /=
        (1 - gamma_.internalField());
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

    UfSmoothed_.internalField() =
        Uf_.internalField()*(1 - gamma_.internalField());

    smoothField(UfSmoothed_);

    UfSmoothed_.internalField() /=
        (1 - gamma_.internalField());

    // evolve Ns steps forward each time when Lammps is called.
    for (label k = 0; k < Ns; k++)
    {
        label nLocal = size(); // Get number of local particles

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
        lammpsEvolveForward(XLocal, VLocal, lmpCpuIdLocal, pDrag_, nstep);

        // update position/velocity of all particles in this cloud.
        // (Harvest XLocal & VLocal)  Lammps --> Cloud
        setPositionVeloCpuId(XLocal, VLocal, lmpCpuIdLocal);

        diffusionRunTime_.cpuTimeIncrement();
        // move particle to the new position
        Cloud<softParticle>::move(td0, mesh_.time().deltaTValue());

        particleMoveTime_ += diffusionRunTime_.cpuTimeIncrement();

        if (particleCount_ != size())
        {
            Pout<< "Warning: enhancedCloud::evolve: "
                << "particle number modified! "
                << "Particle number before: " << particleCount_
                << "Particle number now: " << size()
                << endl;

            particleCount_ = size();
        }

        // make sure all particles are in cell.
        assertParticleInCell();

        // update Uri (relative particle velocities) here.
        updateParticleUr();

        // change Eulerian (mesh-based) alpha field
        particleToEulerianField();

        delete [] XLocal;
        delete [] VLocal;
        delete [] lmpCpuIdLocal;
    }

    Pout<< "After this cycle, "
        << size() << " local particles has been moved. " << endl;
}



//- comments
//- comments
// template <class valueType>
// void enhancedCloud::smoothField(GeometricField <valueType, fvPatchField, volMesh> & sFieldIn)
// {
//   valueType  diffWorkField
//     (
//        IOobject
//          (
//             "tempDiffu",
//             diffusionRunTime_.timeName(),
//             diffusionMesh_,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//        diffusionMesh_,
//        sFieldIn.dimension(),
//        sFieldIn.internalField(),
//        zeroGradientFvPatchScalarField::typeName
//      );
//
//   dimensionedScalar DT("DT", dimensionSet(0, 2, -1, 0, 0), 1.0);
//
//   solve( fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField) );
//
//   sFieldIn.internalField() = diffWorkField.internalField();
// }

void enhancedCloud::smoothField(volScalarField& sFieldIn)
{
    diffusionRunTime_.cpuTimeIncrement();
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

    dimensionedScalar DT("DT", dimensionSet(0, 2, -1, 0, 0), 1.0);

    scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

    Info<< "smoothing " << sFieldIn.name() << endl;

    while (diffusionRunTime_.loop())
    {
        while (simple_.correctNonOrthogonal())
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }
    diffusionRunTime_.setTime(startTime,startIndex);

    sFieldIn.internalField() = diffWorkField.internalField();
    diffusionTimeCount_ += diffusionRunTime_.cpuTimeIncrement();
}


void enhancedCloud::smoothField(volVectorField& sFieldIn)
{
    diffusionRunTime_.cpuTimeIncrement();
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

    dimensionedScalar DT("DT", dimensionSet(0, 2, -1, 0, 0), 1.0);

    scalar startTime = diffusionRunTime_.startTimeIndex();
    label startIndex = diffusionRunTime_.timeIndex();

    Info<< "smoothing " << sFieldIn.name() << endl;

    while (diffusionRunTime_.loop())
    {
        while (simple_.correctNonOrthogonal())
        {
            solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT, diffWorkField));
        }
    }
    diffusionRunTime_.setTime(startTime,startIndex);

    sFieldIn.internalField() = diffWorkField.internalField();
    diffusionTimeCount_ += diffusionRunTime_.cpuTimeIncrement();
}


//- Refresh Ue and Gamma using Gaussian averaging
void enhancedCloud::particleToEulerianField()
{
    gamma_.internalField() *= 0.0;
    Ue_.internalField() *= 0.0;

    // obtain alpha and Ua field
    forAllIter(softParticleCloud, *this, iter)
    {
        softParticle& p = iter();

        label cellI = p.cell();

        // alpha field
        gamma_.internalField()[cellI] +=
            p.Vol();

        Ue_.internalField()[cellI] +=
            p.Vol()*p.U();
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
    smoothField(gamma_);
    smoothField(Ue_);

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
        Utotal2 += Ue_.internalField()[ceI]*mesh_.V()[ceI]*gamma_.internalField()[ceI];
    }

    reduce(Utotal1, sumOp<vector>());
    reduce(Utotal2, sumOp<vector>());

    Info<< "total U solid before: " << Utotal1 << endl;
    Info<< "total U solid after: " << Utotal2 << endl;
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

        if (p.cell()<0) continue;
        dragSumField_.internalField()[p.cell()] += pDrag_[particleI];
    }

    // convert to force on unit cell volume
    dragSumField_.internalField() /= mesh_.V();

    Info<< "Sum up drag on particle in unit cell: "
        << dragSumField_.internalField()
        << endl;
}

} // End namespace Foam


// ************************************************************************* //
