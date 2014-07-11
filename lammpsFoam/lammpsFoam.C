/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Application
    lammpsFoam

Description
    Solver for a system of an incompressible fluid phase with one
    phase dispersed, i.e. particles in a liquid. The dispersed phase
    is modeled with DEM (LAMMPS)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"
#include "Switch.H"
#include "mpi.h"

#include "enhancedCloud.H"
#include "chPressureGrad.H"

//#define DEBUG_HYBRID
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    if (! Pstream::parRun()) MPI_Init(&argc,&argv);

    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"
    scalar t0 = runTime.elapsedCpuTime();
    #include "createParticles.H"
    #include "initContinuityErrs.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalarList splitTime(5,0.0);

    Info<< "\nStarting time loop\n" << endl;
    #include "liftDragCoeffs.H"
        
    splitTime[1] += runTime.elapsedCpuTime() - t0;

    while (runTime.run())
    {
        t0 = runTime.elapsedCpuTime();

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISO.H"
        #include "CourantNo.H"

        #include "alphaEqn.H"

        #include "UEqns.H"


        // --- PISO loop
        #include "pEqn.H"

        #include "DDtU.H"
        #include "kEpsilon.H"

        splitTime[0] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();

        // get drag from latest velocity fields and evolve particles.
        #include "moveParticles.H"

        splitTime[1] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();

        #include "liftDragCoeffs.H"
        #include "write.H"

        splitTime[2] += runTime.elapsedCpuTime() - t0;
        t0 = runTime.elapsedCpuTime();

        #include "writeCPUTime.H"
    }

    Info<< "End\n" << endl;

    if (! Pstream::parRun())  MPI_Finalize();
    return(0);
}


// ************************************************************************* //
