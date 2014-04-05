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
    meshForWeighting

Description
    Construct a mesh with vertex neighbour list
    (neighbour is definded as cellCells sharing vertices)

SourceFiles
    meshForWeighting

\*---------------------------------------------------------------------------*/

#ifndef meshForWeighting_H
#define meshForWeighting_H

// #include "fvMesh.H"
#include "fvCFD.H"
#include "labelList.H"
#include <math.h>
#include <algorithm>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshForWeighting Declaration
\*---------------------------------------------------------------------------*/

class meshForWeighting
{
    // Privete member data

        // - References to the database and meshes
        const fvMesh& mesh_;

        // - Data (vetex cellCells pointer)
        mutable labelListList* vccPtr_;

        // - Calculated according to mesh size
        scalar bandWidth_;

        // - boxCar mask applied
        bool boxCar_;

        // - Bandwidth-CellSize Rario
        // - Should be read from dictionary
        scalar bwCellSizeRatio_;


    // Private memeber functions

        // - Calculating band width
        void calcBandWidth();

        // - Calculating vertex neighbours
        void calcVertexCellCells() const;


protected:

    // - This is a Gaussian-like kernel function
    inline scalar kernel(scalar d)
    {
        return (abs(d) < bandWidth_) ?
            pow(pow(1-pow((d/bandWidth_), 2), 2), 2) : 0.0;
    }

    //- Access vertex neighbors
    const labelListList& vertexCellCells() const;


public:

    //- Construct from mesh (the only way)
    meshForWeighting(const fvMesh& mesh, scalar bwRatio=2.0);

    // Destructor
    ~meshForWeighting();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
