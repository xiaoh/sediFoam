/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::dragModel

Description

SourceFiles
    dragModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
// #include "scalarList.H"
#include "scalarField.H"
#include "dimensionedScalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
{
protected:

    // Protected data

        const dictionary& cloudDict_;

        const IOdictionary& transDict_;

        const scalarField& alpha_;

        const scalarField& pd_;

        scalar nuf_;

        scalar rhof_;

public:

    //- Runtime type information
    TypeName("dragModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            dictionary,
            (
                const dictionary& cloudDict,
                const IOdictionary& transDict,
                const scalarField& alpha,
                const scalarField& pd
            ),
            (
                cloudDict,
                transDict,
                alpha,
                pd
            )
        );


    // Constructors

        dragModel
        (
            const dictionary& cloudDict,
            const IOdictionary& transDict,
            const scalarField& alpha,
            const scalarField& pd
        );


    //- Destructor
    virtual ~dragModel();


    // Selectors

        static autoPtr<dragModel> New
        (
            const dictionary& cloudDict,
            const IOdictionary& transDict,
            const scalarField& alpha,
            const scalarField& pd
        );


    // Member Functions

        //- the dragfunction K used in the momentum eq.
        //    ddt(alpha*rhoa*Ua) + ... = ... alpha*beta*K*(Ua-Ub)
        //    ddt(beta*rhob*Ub) + ...  = ... alpha*beta*K*(Ub-Ua)
        // ********************************** NB ! *****************************
        // for numerical reasons alpha and beta has been
        // extracted from the dragFunction K,
        // so you MUST divide K by alpha*beta when implemnting the drag function
        // ********************************** NB ! *****************************
    virtual tmp<scalarField> Jd(const scalarField& Ur) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
