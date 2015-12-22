/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    label patchID = mesh.boundaryMesh().findPatchID(surfaceName);

    OFstream faceFile1("faceList");

    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

    Info<< "start point is:" << startPointList << endl;
    Info<< "end point is:" << endPointList << endl;
    Info<< "surface name is: " << surfaceName << endl;

    int nPatch = 0;
    forAll(mesh.boundaryMesh()[patchID],faceI)
    {
        for(int i=0;i<2;i++)
        {
            if((cPatch.faceCentres()[faceI].x() - startPointList[i].x())* \
               (cPatch.faceCentres()[faceI].x() - endPointList[i].x()) <= 0 && \
               (cPatch.faceCentres()[faceI].y() - startPointList[i].y())* \
               (cPatch.faceCentres()[faceI].y() - endPointList[i].y()) <= 0 && \
               (cPatch.faceCentres()[faceI].z() - startPointList[i].z())* \
               (cPatch.faceCentres()[faceI].z() - endPointList[i].z()) <= 0)
            {
                faceFile1<< cPatch.start() + faceI << endl;
                nPatch++;
            }
        }
    }

    Info<< "Find " << nPatch << " faces on the patch: " << surfaceName << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
