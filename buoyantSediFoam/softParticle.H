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

Class
    softParticle

Description
    A deformable particle in molecular dynamics.

SourceFiles:
    softParticleI.H
    softParticle.C
    softParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef softParticle_H
#define softParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Class forward declarations
class softParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class softParticle Declaration
\*---------------------------------------------------------------------------*/

class softParticle
:
    public particle
{
    // Private member data

        //- Diameter
        scalar d_;

        //- Mass
        scalar mass_;

        //- Velocity of the particle obtained from LAMMPS.
        //  REAL velocity of the particle in the current time step.
        vector U_;

        //- Velocity used to move the particle in OpenFOAM.
        //  Used when the particle moves across a patch.
        //  Calculated by averaging velocity between two time steps.
        vector moveU_;

        //- Velocity used to for ensemble.
        vector ensembleU_;

        //- Position of the particle in the previous time step.
        //  Used to predict moveU_.
        vector positionOld_;

        //- Velocity of the particle in the previous time step.
        //  Used to calculate the acceleration and compute added mass force.
        vector UOld_;

        //- Particle ID used in LAMMPS to distinguish different particles.
        //  Used as a global particle ID.
        label tag_;

        //- Particle's last seen Lammps processor ID
        label lmpCpuId_;

        //- Particle type. This variable is used to tell different
        //  type of particles (e.g., large or small, heavy or light).
        label type_;

        //- Density of the particle.
        scalar density_;


    // Private Member Functions

        //    //- Disallow default bitwise copy construct
        //        softParticle(const softParticle&);

        //- Disallow default bitwise assignment
        void operator=(const softParticle&);


public:

    friend class Cloud<softParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<softParticleCloud>
    {

    public:

        // Constructors

            trackingData(softParticleCloud& cloud)
            :
                particle::TrackingData<softParticleCloud>(cloud)
            {}
    };


    // Constructors

        //- Construct from components
        softParticle
        (
		     const polyMesh& mesh,
		     const vector& position,
		     const label celli,
		     const scalar& d,
		     const vector& U,
             const scalar& rhos,
             const label& tag,
             const label& lmpCpuId,
             const label& type
		);

        //- Construct from Istream reading field values if required
        softParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new softParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
	        :
	            mesh_(mesh)
            {}

            autoPtr<softParticle> operator()(Istream& is) const
            {
	            return autoPtr<softParticle>(new softParticle(mesh_, is, true));
            }
        };


    // Destructor
    ~softParticle();


    // Member Functions

        // Calculating

            //- Calculate important properties
            void calculateDerived();


        // Tracking

            //- Move
            bool move(trackingData&, const scalar);


        // Access

            //- Return diameter
            scalar& d()
            {
                return d_;
            }

            //- Return mass
            scalar& m()
            {
                return mass_;
            }

            //- Return velocity
            vector& U()
            {
                return U_;
            }

            //- Return moving velocity
            vector& moveU()
            {
                return moveU_;
            }

            //- Return ensemble velocity
            vector& ensembleU()
            {
                return ensembleU_;
            }

            //- Return position in the previous time step
            vector& positionOld()
            {
                return positionOld_;
            }

            //- Return velocity in the previous time step
            vector& UOld()
            {
                return UOld_;
            }

            //- Return global ID
            label& ptag()
            {
                return tag_;
            }

            //- Return global ID
            label& pLmpCpuId()
            {
                return lmpCpuId_;
            }

            //- Return particle type
            label& ptype()
            {
                return type_;
            }

            //- Return volume
            scalar Vol()
            {
                return constant::mathematical::pi*d_*d_*d_/6.0;
            };


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);



    // I-O

        static void readFields(Cloud<softParticle>& c);

        static void writeFields(const Cloud<softParticle>& c);

        static void writeFields(const Cloud<softParticle>& c, const label np);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const softParticle&);
};


template<>
inline bool contiguous<softParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
