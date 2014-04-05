/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-type case

Description
    Test if data sharing OpenFOAM and Lammps can be achieved by casting pointer 
    types without data copying. 
    (foam => lammps: drag force; lammps => foam: particle information X, V etc.) 

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "vector.H"
#include "vectorList.H"
#include "contiguous.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"


    label n = 5;

    vectorList vecs(n);

    for (label i = 0; i < n; i++)
    {
        vecs[i] = vector(i*3, i*3 + 1, i*3 + 2);
    }


    Info << "vecs = " << vecs << endl;

    if (contiguous<vector>())
    {
        Info << "No need to copy memory. Cast pointer instead!" << endl;

        double * array = reinterpret_cast <double *> (&(vecs.first()));

        Info << "array interpretation: " << endl;
        for (label j = 0; j < n*3; j++)
        {
            Info <<  array[j] << endl;
        }
    }
    else
    {
        Info << "You must copy the list instead!" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
