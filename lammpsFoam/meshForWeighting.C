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

\*----------------------------------------------------------------------------*/

#include "meshForWeighting.H"

#include <algorithm>
#include <vector>
namespace Foam
{

using namespace std;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// #define debugVN
void meshForWeighting::calcVertexCellCells() const
{
    if (vccPtr_)
    {
        FatalErrorIn("meshForWeighting::calcVertexCellCells() const")
            << "vertexCellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        // create storage for neighbour list
        vccPtr_ = new labelListList(mesh_.nCells());
        labelListList& cellCellAddr = *vccPtr_;

        forAll(cellCellAddr, cellI)
        {
            std::vector<label> neighbours;

            // push all the neighbours (non-unique for now)
            forAll(mesh_.cellPoints()[cellI], pointI)
            {
                // pID: node ID of each vertex of this cellI
                label pID = mesh_.cellPoints()[cellI][pointI];

                forAll(mesh_.pointCells()[pID], neiI)
                {
                    // cID: cell ID sharing this vertex (pID) with
                    // this cell (cellI)
                    label cID = mesh_.pointCells()[pID][neiI];

                    // don't count itself as a "neighbour".
                    if (cID != cellI)
                    {
                        neighbours.push_back(cID);
                    }
                }
            }

            // sort the list and make it unique
            std::sort(neighbours.begin(), neighbours.end());
            std::vector<label>::iterator endUnique =
                std::unique(neighbours.begin(), neighbours.end());

            neighbours.erase(endUnique, neighbours.end());

            // return the number of neighbours of cellI
            std::vector<label>::size_type ncr = neighbours.size();

            label nci = ncr;
            cellCellAddr[cellI].setSize(nci);

            // push the cell label to cellCellAddr[cellI]
            label i = 0;
            for
            (
                std::vector<label>::const_iterator iter = neighbours.begin();
                iter != neighbours.end();
                ++iter
            )
            {
                cellCellAddr[cellI][i++] = *iter;
            }
        }
    }

#ifdef debugVN
    Info<< "Vertex Neighbours: " << *vccPtr_ << endl;
#endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
meshForWeighting::meshForWeighting
(
    const fvMesh& mesh,
    scalar bwRatio
)
:
    mesh_(mesh),
    bwCellSizeRatio_(bwRatio)
{
    vccPtr_ = NULL;
    calcVertexCellCells();

    // hard-coded lower threshold of bwRatio
    scalar minRatio(0.2);

    if (bwRatio < minRatio)
        boxCar_ = true;
    else
        boxCar_ = false;

    bwRatio = max(bwRatio, minRatio);

    // if (Pstream::parRun())  bwCellSizeRatio_ = 1.5;
    calcBandWidth();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshForWeighting::~meshForWeighting()
{
    if (vccPtr_)
    {
        delete vccPtr_;
        vccPtr_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Access functions
const labelListList& meshForWeighting::vertexCellCells() const
{
    if (!vccPtr_)
    {
        calcVertexCellCells();
    }

    return *vccPtr_;
}


//- Other Operations

//- Calcuate bandwidth
//  The overall mean distance between cells and their face neighbours
//  Note: Here we are not talking about vertex neighbours.
void meshForWeighting::calcBandWidth()
{
    label nc = 0;
    scalar sumDist = 0.0;

    // calculate the total distance
    forAll(mesh_.cells(), cellI)
    {
        scalar meanDist = 0.0;
        forAll(mesh_.cellCells()[cellI], neiI)
        {
            label neiCID = mesh_.cellCells()[cellI][neiI];
            meanDist += mag(mesh_.C()[cellI] - mesh_.C()[neiCID]);
        }

        meanDist /= scalar(mesh_.cellCells()[cellI].size());
        sumDist += meanDist;
        nc++;
    }

    // calculate the band width
    if (nc > 1)
    {
        bandWidth_ =  bwCellSizeRatio_*sumDist/scalar(nc);

        if (boxCar_)
        {
            Info<< "*** NOTE: box-car mask is applied. \n"
                << "*** Bandwidth printed below is not actual!" << endl;
        }

        Info<< "Bandwidth used in ensembling is: "
            << bandWidth_  << ", "
            << " bandwidth/cell size ratio is: " << bwCellSizeRatio_
            <<  endl;
    }
    else
    {
        FatalErrorIn("meshForWeighting::calcBandWidth ")
            << "Cell number is zero\n"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
