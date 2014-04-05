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

\*---------------------------------------------------------------------------*/

#include "SyamlalOBrien.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SyamlalOBrien, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        SyamlalOBrien,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SyamlalOBrien::SyamlalOBrien
(
    const dictionary& cloudDict,
    const IOdictionary& transDict,
    const scalarField& alpha,
    const scalarField& pd
)
:
    dragModel
    (
        cloudDict,
        transDict,
        alpha,
        pd
    )
{
    dimensionedScalar Dnuf_(transDict_.lookup("nub"));
    dimensionedScalar Drhof_(transDict_.lookup("rhob"));

    nuf_ = Dnuf_.value();
    rhof_ = Drhof_.value();

    Info<< "Fluid properties dictionary look-up from Drag model: " << endl
        << "Kinematic Viscosity: " << Dnuf_ << "; " << "Density: "
        << Drhof_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SyamlalOBrien::~SyamlalOBrien()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//  #define DEBUG_JD

Foam::tmp<Foam::scalarField> Foam::SyamlalOBrien::Jd
(
    const scalarField& Ur
) const
{

    if
    (
        Ur.size() != alpha_.size()
     || Ur.size() != pd_.size()
    )
    {
        FatalErrorIn("SyamlalOBrien::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;
    }

    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));
    scalarField Ai = pow(beta, 4.14);
    scalarField Bi = 0.8*pow(beta, 1.28);

    forAll(Ur, pI)
    {
        if (beta[pI] > 0.85)
        {
            Bi[pI] = pow(beta[pI], 2.65);
        }
    }

    scalarField Re = max(Ur*pd_/nuf_, scalar(1.0e-3));

    scalarField Vr =
        0.5
       *(
            Ai - 0.06*Re + sqrt(sqr(0.06*Re)
          + 0.12*Re*(2.0*Bi - Ai) + sqr(Ai))
        );

    scalarField Cds = sqr(0.63 + 4.8*sqrt(Vr/Re));

#ifdef DEBUG_JD
    Info<< " ==Report====> " << "SyamlalOBrien::Jd()  " << endl;
    Info<< "Lagrangian A: " << Ai
        << "Lagrangian B:" << Bi
        << "Lagrangian rho: " << rhof_
        << "Lagrangian nu: " << nuf_
        << "Lagrangian pd: " << pd_
        << "Lagrangian alpha: " << alpha_
        << "Lagrangian Re: " << Re
        << "Lagrangian Ur:" << Ur
        << "Lagrangian Vr: " << Vr
        << "Lagrangian Cds: " << Cds
        << "Lagrangian Jd: " << 0.75*Cds*rhof_*Ur/(pd_*sqr(Vr)) << endl;
#endif

      return 0.75*Cds*rhof_*Ur/(pd_*sqr(Vr));
}

// ************************************************************************* //
