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
    Foam::NoCorrection

Description
    Syamlal, M., Rogers, W. and O'Brien, T. J. (1993) MFIX documentation,
    Theory Guide. Technical Note DOE/METC-94/1004. Morgantown, West Virginia,
    USA.

SourceFiles
    NoCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef NoCorrection_H
#define NoCorrection_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NoCorrection Declaration
\*---------------------------------------------------------------------------*/

class NoCorrection
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("NoCorrection");


    // Constructors

        //- Construct from components
        NoCorrection
        (
            const dictionary& cloudDict,
            const IOdictionary& transDict,
            const scalarField& alpha,
            const scalarField& pd
        );


    //- Destructor
    virtual ~NoCorrection();


    // Member Functions

        tmp<scalarField> Jd(const scalarField& Ur) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
