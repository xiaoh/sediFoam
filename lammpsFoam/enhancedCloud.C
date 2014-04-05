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

//- Update weights corresponding to all particles.
void enhancedCloud::updateCellWeights()
{
    if (!weightPtr_)
    {
        // create storage for weight list if not allocated yet
        weightPtr_ = new scalarListList(particleCount_);
    }

    if (!particleCellPtr_)
    {
        // create storage for cell list if not allocated yet
        particleCellPtr_ = new labelListList(particleCount_);
    }

    // if particle number modified
    if (particleCount_ != size())
    {
        Pout<< "Warning: enhancedCloud::updateCellWeights: "
            << "softParticleCloud object modified! "
            << " Particle count before: " << particleCount_
            << " Particle count now: " << size()
            << endl;
        particleCount_ = size();
        reAllocateLists(); // Re-allocate lists
    }

    scalarListList& pCellWeights = *weightPtr_;
    labelListList& pCellLabels = *particleCellPtr_;

    int particleI = 0;
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter, ++particleI
    )
    {
        softParticle& p = pIter();
        label cellP = p.cell();

        if (cellP<0) continue;

        const labelList& neiCellIDs = vertexCellCells()[cellP];
        label neiN = neiCellIDs.size();

        pCellWeights[particleI].setSize(neiN+1);
        pCellLabels[particleI].setSize(neiN+1);

        // update the neighbor cells of this particle
        forAll(neiCellIDs, neiI)
        {
            label ncID = neiCellIDs[neiI];
            pCellLabels[particleI][neiI] = ncID;
        }

        // finally, the cell associated with this particle.
        pCellLabels[particleI][neiN] = cellP;

        // calculate normalized weights of cells for this particle
        calcWeightGroup
        (
            pCellWeights[particleI],
            p,
            cellP,
            particleI
        );
    }
}


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
        if (cellI<0) continue;
        // pAlpha_[particleI] = gamma_[cellI];
        pAlpha_[particleI] = ensembleAlphaTimeFixed_[cellI];
    }
}


//- Update relative Ur list (per fluid step or substep)
//  Update Ur_ and magUr_
//  Called by calcTcField before calculating Omega & Asrc
//  Also called after each substep.
void enhancedCloud::updateParticleUr()
{
    softParticleCloud::iterator pIter = begin();

    Uri_.setSize(particleCount_);
    magUri_.setSize(particleCount_);

    forAll(particleWeights(), particleI)
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
        vector wtUfi(vector::zero);

        /*
            forAll(particleWeights()[particleI], neiI)
            {
                label cellI = particleCells()[particleI][neiI];
                scalar weight = particleWeights()[particleI][neiI];
                #ifdef DEBUG_ENH
                Info<< "cellI: " << cellI
                    << " Ucell: " << Uf_[cellI]
                    << " weight: " << weight
                    << endl;
                #endif
                wtUfi += weight*Uf_[cellI];
            }
        */

        wtUfi = Uf_[p.cell()];
        Uri_[particleI] = wtUfi - p.U();
        magUri_[particleI] = mag(Uri_[particleI]);

        ++pIter;

#ifdef DEBUG_ENH
        Info<< "particle U: " << p.U()
            << " weight Ufi: " <<  wtUfi
            <<" Uri: " << Uri_[particleI] << endl;
#endif
    }

    if (pIter != end())
    {
        FatalErrorIn("enhancedCloud::updateParticleUr()")
            << "Inconsistent particle correspondance."
            << abort(FatalError);
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

    // scan particle list to average and get Omega & Asrc field
    softParticleCloud::iterator pIter = begin();
    forAll(particleWeights(), particleI)
    {
        softParticle& p = pIter();

        // update Omega field (fluid density omitted)
        label cellI = p.cell();

        if (cellI < 0) continue;

        scalar omg = p.Vol()*Jd_[particleI]/(mesh_.V()[cellI]);
        Omega_.internalField()[cellI] += omg;
        Asrc_.internalField()[cellI] += omg*p.ensembleU();
        ++pIter;
        /*
            forAll(particleWeights()[particleI], neiI)
            {
                label cellI = particleCells()[particleI][neiI];
                scalar weight = particleWeights()[particleI][neiI];

                // Update Omega field (density omitted)
                scalar omg = weight*
                    p.Vol()*Jd_[particleI]/(mesh_.V()[p.cell()]);

                // This is still debating ...
                Omega_.internalField()[cellI] += omg;
                Asrc_.internalField()[cellI] += omg * ensPU_[particleI];
            }
        */
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
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
    scalar bwDxRatio
)
:
    softParticleCloud(vpi, U, p, Ue, nu, alpha, cloudDict),
    meshForWeighting(U.mesh(), bwDxRatio),
    mesh_(U.mesh()),
    Uf_(Uf),
    Omega_
    (
        IOobject
        (
            "Omega",
            runTime().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
           "zero",
           dimensionSet(1, -2, -2, 0, 0),
           vector::zero
        )
    )
{
    weightPtr_ = NULL;
    particleCellPtr_ = NULL;

    drag_ = Foam::dragModel::New(cloudDict, transDict, pAlpha_, pDia_);

    particleCount_ = size();

    // a threshhold to turn on box-car option
    // don't change the hard-coded threhold value!
    if (bwDxRatio < 0.2)
        boxCar_ = true;
    else
        boxCar_ = false;

    // initial quantities for dragModel
    pDia_.setSize(particleCount_);
    pAlpha_.setSize(particleCount_);
    Uri_.setSize(particleCount_);
    magUri_.setSize(particleCount_);

    // initialise drag force on each particle
    pDrag_ = vectorList(particleCount_);

    // initialise ensemble velocity on each particle
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter
    )
    {
        softParticle& p = pIter();
        p.ensembleU() = p.U();
    }

    // initialize cell weights
    updateCellWeights();

    // initialize alpha field and Ue
    particleToEulerianField();

    // setup ensAlpha_
    ensembleAlpha_ = gamma_.internalField();
    ensembleAlphaTimeFixed_ = ensembleAlpha_;

    setupParticleDia();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

enhancedCloud::~enhancedCloud()
{
    if (weightPtr_)
    {
        delete weightPtr_;
        weightPtr_ = 0;
    }

    if (particleCellPtr_)
    {
        delete particleCellPtr_;
        particleCellPtr_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void enhancedCloud::evolve()
{
    softParticle::trackingData td0(*this);

    label Ns = subCycles_;

    initEnsemble();

    // evolve Ns steps forward each time when Lammps is called.
    for (label k = 0; k < Ns; k++)
    {
        label nLocal = size(); // Get number of local particles

        // the squence of the data on XLocal and VLocal is the same as
        // sequence of local particle index.
        vector* XLocal = new vector[nLocal];
        vector* VLocal = new vector[nLocal];

        int nstep = subSteps_;

        // here, Ur (relative velocity is called)
        updateParticleUr();

        updateDragOnParticles();

        // XLocal & VLocal are work spaces for "lammpsEvolveForward"
        // newly obtianed values are put there
        lammpsEvolveForward(XLocal, VLocal, pDrag_, nstep);

        // update position/velocity of all particles in this cloud.
        // (Harvest XLocal & VLocal)  Lammps --> Cloud
        setPositionVelo(XLocal, VLocal);

        // move particle to the new position
        Cloud<softParticle>::move(td0, mesh_.time().deltaTValue());
        // setPositionCell();

        updateCellWeights();

        // make sure all particles are in cell.
        // assertParticleInCell();

        // update Uri (relative particle velocities) here.
        updateParticleUr();

        // change Eulerian (mesh-based) alpha field
        particleToEulerianField();

        accumulateEnsemble();

        delete [] XLocal;
        delete [] VLocal;
    }

    // compute ensembled value after n iterations
    computeEnsemble();

    Info<< "After this cycle, "
        << size() << " local particles has been moved. " << endl;
}


//- Initialize ensemble average of alpha and particle velocity
void  enhancedCloud::initEnsemble()
{
    // set accumulation count to zero.
    ensembleCount_ = 0;
    ensembleAlpha_ *= 0.0;

    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter
    )
    {
        softParticle& p = pIter();
        p.ensembleU() = vector::zero;
    }
}


//- Accumulate ensemble average of alpha and particle velocity.
//  Warning: using ensemblePU_ and ensembleAlpha_ during sub-stepping is error!
void  enhancedCloud::accumulateEnsemble()
{
    ++ensembleCount_;

    // accumulate Alpha for ensemble
    ensembleAlpha_ += gamma_.internalField();

    // accumulate particle velocity for ensemble
    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter
    )
    {
        softParticle& p = pIter();
        p.ensembleU() += p.U();
    }
}


//- Compute ensemble average of alpha and particle velocity
void  enhancedCloud::computeEnsemble()
{
    ensembleAlpha_ /= ensembleCount_;
    ensembleAlphaTimeFixed_ = ensembleAlpha_;

    for
    (
        softParticleCloud::iterator pIter = softParticleCloud::begin();
        pIter != softParticleCloud::end();
        ++pIter
    )
    {
        softParticle& p = pIter();
        p.ensembleU() /= ensembleCount_;
    }
}


//- Calculate the weight for a group. Should be inline.
//  Called only by calcCellWeights() for each particle.
inline void enhancedCloud::calcWeightGroup
(
    scalarList& weight,
    softParticle& p,
    label& cellI,
    int& i
)
{
    const labelList& neiCellIDs = vertexCellCells()[cellI];
    scalar dist = 0.0;
    scalar totalWt = 0.0;
    label neiN = neiCellIDs.size();

    // short circuit if bwDxRatio is smaller than a threshold.
    if (boxCar_)
    {
        weight = 0.0;
        weight[neiN] = 1.0;
        totalWt = 1.0;
        return;
    }

    forAll(neiCellIDs, neiI)
    {
        label ncID = neiCellIDs[neiI];
        dist = mag(p.position()- mesh_.C()[ncID]);
        weight[neiI] = kernel(dist);
        totalWt += weight[neiI];
    }

    // now, the cell this particle belongs to
    dist = mag(p.position()-mesh_.C()[cellI]);
    weight[neiN] = kernel(dist);
    totalWt += weight[neiN];

    // bandwidth too small for this cell.
    if (totalWt < 1e-4)
    {
        Pout<< "Warning: Bandwidth too small for this cell." << endl
            << " cell ID: " << cellI
            << " Particle ID: " << i
            << " neighbors ID: " << neiCellIDs
            << " particle x: "  << p.position()
            << " associated cell center: " << mesh_.C()[cellI]
            << " Distance: " << dist
            << " weight: " << weight
            << endl;

        // Box-car kernel function of improper band width specified.
        weight = 0.0; weight[neiN] = 1.0; totalWt = 1.0;
    }
    else
    {
        // normalize all weight by total weight
        forAll(neiCellIDs, neiI)
        {
            weight[neiI] /= totalWt;
        }

        weight[neiN] /= totalWt;
    }
}


//- Refresh Ue, and Gamma using Gaussian averaging
void enhancedCloud::particleToEulerianField()
{
    gamma_.internalField() *= 0.0;
    Ue_.internalField() *= 0.0;

    softParticleCloud::iterator pIter = begin();
    forAll(particleWeights(), particleI)
    {
        softParticle& p = pIter();

        if (p.cell() < 0) { ++pIter; continue;}

        forAll(particleWeights()[particleI], neiI)
        {
            label cellI = particleCells()[particleI][neiI];
            scalar weight = particleWeights()[particleI][neiI];

            // alpha field
            gamma_.internalField()[cellI] +=
                p.Vol()*weight;

            // particle velocity field
            Ue_.internalField()[cellI] +=
                p.Vol()*p.U()*weight;

            // should be mass-averaged velocity.
            // but the density data have not been obtained from Lammps yet.
            // use volume average for now.
        }

        ++pIter;
    }

    gamma_.internalField() /= mesh_.V();

    forAll(Ue_.internalField(), ceI)
    {
        if (gamma_.internalField()[ceI] > 1e-10)
        {
            Ue_.internalField()[ceI] /=
                (mesh_.V()[ceI]*gamma_.internalField()[ceI]);
        }
    }

    if (pIter != end())
    {
        FatalErrorIn("enhancedCloud::particleToEulerianField")
            << "Inconsistent particle correspondance."
            << abort(FatalError);
    }
}


const labelListList& enhancedCloud::particleCells()
{
    if (!particleCellPtr_)
    {
        updateCellWeights();
    }

    return *particleCellPtr_;
}


const scalarListList& enhancedCloud::particleWeights()
{
    if (!weightPtr_)
    {
        updateCellWeights();
    }

    return *weightPtr_;
}


//- Reallocate listLists if particle number has changed.
void enhancedCloud::reAllocateLists() const
{
    if (weightPtr_)
    {
        delete weightPtr_;
        weightPtr_ = new scalarListList(particleCount_);
    }

    if (particleCellPtr_)
    {
        delete particleCellPtr_;
        particleCellPtr_ = new labelListList(particleCount_);
    }
}


//- Check if all particles in fluid cells
//  It is in this class because it used to depend on neighbour info.
//  Not any more though.
void enhancedCloud::assertParticleInCell()
{
    softParticleCloud::iterator pIter = begin();

    forAll(particleCells(), particleI)
    {
        softParticle& p = pIter();
        point& pos = p.position();

        bool found = false;

        label lCellI = p.cell();

        if (lCellI >= 0) found = true;

        if (!found)
        {
            WarningIn("enhancedCloud::assertParticleInCell()")
                << "Particle outside fluid domain!" << endl
                << " Debug Info: "
                << " Particle ID: " << particleI
                << " Particle Position: " << pos
                << " Associated cell ID: " << p.cell()
                << abort(FatalError);
        }

        ++pIter;
    }

} // End function


void enhancedCloud::weightInfo()
{
    Info<< "particleCells: " << particleCells() << endl;
    Info<< "particleWeights: " << particleWeights() << endl;
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
