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

#include "ErgunWenYu.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ErgunWenYu, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        ErgunWenYu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ErgunWenYu::ErgunWenYu
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

Foam::ErgunWenYu::~ErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// #define DEBUG_JD3

Foam::tmp<Foam::scalarField> Foam::ErgunWenYu::Jd
(
    const scalarField& Ur
) const
{

    if
    (
        Ur.size() != alpha_.size()
     || Ur.size() != pd_.size()
    )
        FatalErrorIn("ErgunWenYu::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;

    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));
    scalarField bp = pow(beta, -2.65);
    scalarField Re = max(beta*Ur*pd_/nuf_, scalar(1.0e-3));
    scalarField Cds = 24.0*(1.0 + 0.15*pow(Re, 0.687))/Re;

    forAll(Re, particleI)
    {
        if (Re[particleI] > 1000.0)
        {
            Cds[particleI] = 0.44;
        }
    }

    // Wen and Yu (1966)
    tmp<scalarField> tKWenYu = 0.75*Cds*rhof_*Ur*bp/pd_;
    scalarField& KWenYu = tKWenYu();

    // Ergun
    forAll(beta, particleJ)
    {
        if (beta[particleJ] <= 0.8)
        {
            KWenYu[particleJ] =
                150.0*alpha_[particleJ]*nuf_*rhof_
               /sqr(beta[particleJ]*pd_[particleJ])
              + 1.75*rhof_*Ur[particleJ]
               /(beta[particleJ]*pd_[particleJ]);
        }
    }

#ifdef DEBUG_JD3
    Info<< " ==Report====> " << "ErgunWenYu::Jd()  " << endl;
    Info<< " alpha: " << alpha_
        << " Re: " << Re
        << " Ur:" << Ur
        << " Cds: " << Cds
        << " Jd: " << tKWenYu << endl;
#endif

    return tKWenYu;

}

// ************************************************************************* //
