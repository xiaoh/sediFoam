/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

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

Class
    enhancedCloud

Description
    softParticleCloud with averaging features
    and dragModel

SourceFiles
    enhancedCloud.C

\*---------------------------------------------------------------------------*/

#ifndef enhancedCloud_H
#define enhancedCloud_H

#include "softParticle.H"
#include "softParticleCloud.H"
#include "dragModel.H"
#include "vectorList.H"
#include "fvPatchField.H"
#include "volMesh.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class template Declaration
\*---------------------------------------------------------------------------*/

class enhancedCloud
:
    public softParticleCloud
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Fluid velocity
        const volVectorField& Uf_;

        //- Smoothed fluid velocity
        volVectorField UfSmoothed_;

        //- Number of particles as remembered by weight operations
        label particleCount_;

        //- List of alpha corresponding to each particle
        scalarField pAlpha_;

        //- List of particle diameter
        scalarField pDia_;

        //- List of particle drag force
        vectorList pDrag_;

        //- Jd coefficient
        scalarField Jd_;

        //- Drag model
        autoPtr<Foam::dragModel> drag_;

        //- List of relative velocity (mag) of each particle
        vectorField Uri_;
        scalarField magUri_;

        //- Eulerian Fields for Tc
        volScalarField Omega_;
        volVectorField Asrc_;

        //- Fluid density
        // scalar rhob_;

        Time diffusionRunTime_;
        fvMesh diffusionMesh_;
        simpleControl simple_;

        scalarList diffusionTimeCount_;
        scalar particleMoveTime_;

        label UfSmoothFlag_;

    // Private Member Functions

        //- Use diffusion equation to smooth ensembled fields.
        // template <class valueType>
        // void smoothField(GeometricField <valueType, fvPatchField, volMesh> & fieldIn);

        void smoothField(volScalarField& );
        void smoothField(volVectorField& );

        //- Setup particle diameter
        void setupParticleDia();

        //- Update particle alpha list (per fluid step)
        void updateParticleAlpha();

        //- Update relative Ur list
        //- (per fluid step or substep)
        //- Update Ur_ and magUr_
        void updateParticleUr();

        //- Compute averge quantities
        //- particle quantities --> Eulerian quantities
        void particleToEulerianField();

        //- drag on each particle
        //- update pDrag_
        void  updateDragOnParticles();

        // Assert all particles in cells;
        void assertParticleInCell();


public:

    // Constructors

        //- Construct null
        enhancedCloud
        (
             const volVectorField& U,
             const volScalarField& p,
             volVectorField& Ue,
             const volVectorField& Uf,
             dimensionedScalar nu,
             volScalarField& alpha,
             IOdictionary& cloudDict,
             IOdictionary& transDict,
             scalar diffusionBandWidth,
             scalar diffusionSteps
        );


    // Destructor
    ~enhancedCloud();


    // Member Functions

        //- Use averaging method to get Omega and A field
        void calcTcFields();

        //- Evolve function
        void evolve();


        // Access

            //- Return reference of Omega field
            const volScalarField& Omega() const
            {
                return Omega_;
            }

            //- Return reference of A field
            const volVectorField& Asrc() const
            {
                return Asrc_;
            }

            //- Return number of particles in the cloud
            label particleCount() const
            {
                return particleCount_;
            }

            //- Return total diffusion time
            const scalarList& diffusionTimeCount() const
            {
                return diffusionTimeCount_;
            }

            //- Return total diffusion time
            const scalar& particleMoveTime() const
            {
                return particleMoveTime_;
            }

            //- Print particle drag sum-up field
            //--- Not parallel yet.
            void dragInfo();

            //- Print average info of the particles
            void averageInfo();

}; // End of class "enhancedCloud"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
