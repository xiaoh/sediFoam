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
    softParticleCloud

Description
    A Cloud of soft particles

SourceFiles
    softParticleCloudI.H
    softParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef softParticleCloud_H
#define softParticleCloud_H

#include "fvCFD.H"
#include "Cloud.H"
#include "error.H"
#include "Pstream.H"
#include "vectorList.H"
#include "tensorList.H"

#include "LammpsCollection.H"
#include "softParticle.H"
#include "interpolation.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

using  namespace LAMMPS_NS;
using  namespace std;

/*---------------------------------------------------------------------------*\
                       Class softParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class softParticleCloud
:
    public Cloud<softParticle>
{
    // Private data

        //- LAMMPS
        LAMMPS* lmp_;

        // Temporarily holder for particle properties in LAMMPS

            //- Particle position
            double* xArray_;

            //- Particle velocity
            double* vArray_;

            //- Particle diameter
            double* dArray_;

            //- Particle drag force
            double* fArray_;

            //- Particle density
            double* rhoArray_;

            //- Particle ID in LAMMPS
            int* tagArray_;

            //- Processor ID of the particles in LAMMPS
            int* lmpCpuIdArray_;

            //- Particle type
            int* typeArray_;

        //- Number of particles of all processors
        int nGlobal_;

        //- Kinematic viscosity of fluid
        dimensionedScalar nu_;

        // References to the database and meshes

            const Time& runTime_;

            const fvMesh& mesh_;

        //- Cycles in LAMMPS per time step
        //  Relation:  solidStepsPerDt_ = nsubsteps * subCycles + nExtra
        scalar solidStepsPerDt_;

        //- Cpu time spent on different parts
        scalarList cpuTimeSplit_;


    // Private Member Functions

        // Lammps related functions

            //- Initialization of LAMMPS
            void initLammps();

            //- Destructor of LAMMPS
            void finishLammps();

            //- Construct OpenFOAM particles
            void initConstructParticles
            (
                int Np,
                double* x,
                double* v,
                double* d,
                double* rho,
                int* tag,
                int* lmpCpuId,
                int* type
            );

        //- Disallow default bitwise copy construct
        softParticleCloud(const softParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const softParticleCloud&);

        //- Adjust the time step for LAMMPS
        void adjustLampTimestep();

        //- Add OpenFOAM particles
        void addNewParticles();

        //- Add and delete particles
        void addAndDeleteParticle();

        //- find the CFD cell ID in which the particles are added
        void findAddParticleCells();


protected:

    // Protected data

        //- Number of sub-cycles
        scalar subCycles_;

        //- Number of extra steps
        label nXtra_;

        //- Number of steps in a sub-cycle
        label subSteps_;

        //- Number of particles in this processor
        label nLocal_;

        // Reference to the Eulerian (mesh-based) field of particle velocity
        // Obtained by ensemble average from LAMMPS
        volVectorField& Ue_;

        IOdictionary cloudProperties_;

        //- Gravity
        vector gravity_;

        //- If add any particles in the simulation
        scalar addParticleOption_;
        tensor addParticleBox_;
        vector addParticleBoxEccentricity_;
        scalar addParticleTimeStep_;
        scalar timeToAddParticle_;
        labelList addParticleLocalList_;
        labelList addParticleCellID_;
        vector addParticleInfo_;
        scalar randomPerturb_;
        tensorList lmpLocalBoxList_;
        label maxTag_;

        scalar deleteBeforeAddFlag_;
        labelList deleteBeforeAddList_;
        labelList deleteParticleList_;
        tensor inletBox_;
        tensor clearInitialBox_;

        label totalAdd_;
        label totalDelete_;
        label totalDeleteBeforeAdd_;


        //- If delete any particles in the simulation
        scalar deleteParticleOption_;
        tensor deleteParticleBox_;

        // References to the physical field

            //- Velocity
            const volVectorField& U_;

            //- Pressure of fluid
            const volScalarField& pf_;

            //- Solid volume fraction of the cell
            volScalarField& gamma_;


    // Protected Member Functions

        //- Simulate the particles motion in LAMMPS
        void lammpsEvolveForward
        (
            vector* XLocal,
            vector* VLocal,
            int* lmpCpuIdLocal,
            vectorList FLocal,
            int nstep
        );

        //- Set particle positions and velocities of the cloud
        //  using data from LAMMPS
        void setPositionVeloCpuId
        (
            vector* XLocal,
            vector* VLocal,
            int* lmpCpuIdLocal
        );

        //- Set the particle cell index after the particles
        //  move across the processor boundary
        void setPositionCell();

        //- Determine if a point is in the region
        bool pointInRegion(vector& point, tensor& box);

        //- Determine if a point is in the box
        bool pointInBox(vector& point, tensor& box);


        //- The transpose data only for the neighbouring processors
        Switch transposeNbrOnly_;

        //- Transpose the scalar of different processors
        template<class DataType>
        void transposeAmongProcs
        (
            List<DataType>& toEveryone,
            List<DataType>& fromEveryone
        );

        //- Convert a long list into list of lists
        template <class DataType>
        void assembleList
        (
            DataType& inList,
            List<DataType>& outList,
            labelList indexList,
            labelList cpuList
        );

        //- Convert a list of lists to a longer list
        template <class DataType>
        void flattenList
        (
            List<DataType>& inList,
            DataType& outList
        );

public:

    // Constructors

        //- Construct from components
        softParticleCloud
        (
			const volVectorField& U,
			const volScalarField& p,
            volVectorField& Ue,
            dimensionedScalar nu,
            volScalarField& alpha,
            IOdictionary& cloudProperties
		);


    // Destructor
    ~softParticleCloud();


    // Member Functions

        // Access

            //- Return run time
            const Time &runTime() const
            {
                return runTime_;
            }

            //- Return number of sub-cycles
            const scalar &subCycles() const
            {
                return subCycles_;
            }

            //- Return number of local particles
            const label &nLocal() const
            {
                return nLocal_;
            }

            //- Return the size of the particle list
            label size() const
            {
                return IDLList<softParticle>::size();
            };

            //- Return the first iterator of the particle list
            iterator begin()
            {
                return IDLList<softParticle>::begin();
            };

            //- Return the last iterator of the particle list
            iterator end()
            {
                return IDLList<softParticle>::end();
            };

            //- Cpu time spent on different parts
            scalarList cpuTimeSplit()
            {
                return cpuTimeSplit_;
            };

        // I-O

            virtual void writeFields() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
