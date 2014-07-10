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
  chPressureGrad

  Description
  A class for applying pressure gradient for a channel flow to keep
  constant mass flux.

  SourceFiles
  chPressureGrad.H
  chPressureGrad.C
  \*---------------------------------------------------------------------------*/


#ifndef   chPressureGrad_H
#define   chPressureGrad_H

#include "Switch.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class chPressureGrad Declaration
  \*---------------------------------------------------------------------------*/

  class chPressureGrad
  {
    // Private data

  private:
    //- static variables    
    static bool globalParamSet_;
    static vector flowDirection_;
    static Switch chFlowMode_;
    static dimensionedVector Ubar_;
    static dimensionedVector gradPbar_;
    static dimensionedScalar magUbar_;
    static word specifiedQuantity_;

    volVectorField& U_;
    const volScalarField& alpha_;
    const word name_;
    const word solverName_;

    dimensionedScalar value_;
    IOdictionary gradPDict_;
 
    // Private Member Functions

  public:


    // Constructors
    
    //- Construct from components
    chPressureGrad
    (
        volVectorField& U,
        volScalarField& alpha,
        word name, word solverName,
        word dictDir="uniform"
    );


    // Destructor
    ~chPressureGrad(){};

    // Member Functions (public)
    void adjust(const volScalarField& rUA);

    static void  initPressureGrad(const IOdictionary& transportProperties);
    
    // Access
    word name() {return name_;}
    dimensionedScalar value() {return value_;}
    Switch chFlowMode() {return chFlowMode_;}
    vector flowDirection() {return flowDirection_;}

    // Check

    // Edit

    // Write

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
