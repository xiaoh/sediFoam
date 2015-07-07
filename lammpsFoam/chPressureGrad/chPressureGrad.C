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

  \*---------------------------------------------------------------------------*/

#include "chPressureGrad.H"
#include <fstream>

namespace Foam {
  bool chPressureGrad::globalParamSet_ = false;
  Switch chPressureGrad::chFlowMode_ = false;
  vector chPressureGrad::flowDirection_ = vector::zero;
  dimensionedVector  chPressureGrad::Ubar_("Ubar", dimVelocity, vector::zero);
  dimensionedVector  chPressureGrad::gradPbar_("gradPbar", dimAcceleration, vector::zero);
  dimensionedVector  chPressureGrad::varyingGradP_("varyingGradP", dimAcceleration, vector::zero);
  dimensionedVector  chPressureGrad::dpdt_("dpdt", dimAcceleration, vector::zero);
  dimensionedScalar  chPressureGrad::magUbar_ = mag(Ubar_);
  dimensionedScalar  chPressureGrad::period_("period", dimTime, 0.0);
  scalarList  chPressureGrad::pressureList_(0);

  word chPressureGrad::varyingType_("none");
  word chPressureGrad::specifiedQuantity_("none");

  void  chPressureGrad::initPressureGrad(const IOdictionary& transportProperties)
  {
    
      if(transportProperties.found("Ubar"))
      {
          if(transportProperties.found("gradPbar"))
          {
              FatalErrorIn("chPressureGrad::initPressureGrad()") 
                  << nl
                  << "You specified both Ubar and gradPbar in trasnportProperties dictionary." << nl
                  << "Set only one of them!"<< endl
                  << abort(FatalError);
          }
          
          //  Read centerline velocity for channel simulations
          Ubar_ = transportProperties.lookup("Ubar");
          magUbar_ = mag(Ubar_);
          flowDirection_ = (Ubar_/magUbar_).value();
          specifiedQuantity_ = word("Ubar");
          chFlowMode_ = true;
          
          Info << "*********** NOTE ***********" << nl
               << "Running in CHANNEL FLOW MODE. " << nl
               << "Constant mass flux enforced!" << nl
               << "Ubar = " << Ubar_.value() << nl 
               << "****************************" << nl
               << endl;
      }
      else if(transportProperties.found("gradPbar"))
      {
          gradPbar_ = transportProperties.lookup("gradPbar");
          specifiedQuantity_ = word("gradPbar");
          dimensionedScalar magGradPbar = mag(gradPbar_);
          magGradPbar.value() += VSMALL;
          flowDirection_ = (gradPbar_ /magGradPbar ).value();
          chFlowMode_ = true;

          Info << "*********** NOTE ***********" << nl
               << "Running in CHANNEL FLOW MODE. " << nl
               << "constant gradP imposed!" << nl
               << "gradP = " << gradPbar_.value() << nl 
               << "****************************" << nl
               << endl;

          if(transportProperties.found("dpdt"))
          {
              dpdt_ = transportProperties.lookup("dpdt");
              Info << "*********** NOTE ***********" << nl
                   << "gradP will increase by dpdt!" << nl
                   << "dpdt = " << dpdt_.value() << nl 
                   << "****************************" << nl
                   << endl;
          }
          
      }
      else if(transportProperties.found("varyingGradP"))
      {
          varyingGradP_ = transportProperties.lookup("varyingGradP");
          period_ = transportProperties.lookup("varyingPeriod");
          Foam::word varyingType = transportProperties.lookup("varyingType");
          varyingType_ = varyingType;
          specifiedQuantity_ = word("varyingGradP");
          dimensionedScalar maxGradP = mag(varyingGradP_);
          maxGradP.value() += VSMALL;
          flowDirection_ = (varyingGradP_ /maxGradP ).value();
          chFlowMode_ = true;

          std::ifstream ifs;
          if (varyingType_ == "input")
          {
              double data;
              ifs.open("in.pressure");
              int i = 1;
              while (ifs >> data)
              {
                  pressureList_.setSize(i);
                  pressureList_[i - 1] = data;
                  i++;
              }
          }


          Info << "*********** NOTE ***********" << nl
               << "Running in CHANNEL FLOW MODE. " << nl
               << "varying gradP imposed!" << nl
               << "varying type " << varyingType_ << nl
               << "maxP = " << maxGradP.value() << nl 
               << "****************************" << nl
               << endl;
      }
      else
      {
        Info << "+++++++++++ NOTE +++++++++++++"<< nl
             << "Running in REGULAR MODE."      << nl
             << "NO PRESSURE GRADIENT APPLIED!" << nl
             << "++++++++++++++++++++++++++++++"
             << endl;
      }
   
    globalParamSet_ = true;
  }

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
  // * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  chPressureGrad::chPressureGrad
  (
   volVectorField& U,
   volScalarField& alpha,
   word name,
   word solverName,
   word dictDir
   )
  :
    U_(U), 
    alpha_(alpha), 
    name_(name),
    solverName_(solverName),
    value_(name,  dimensionSet(0, 1, -2, 0, 0), 0.0),
    gradPDict_
    (
     IOobject
     (
      name,
      U.time().timeName(),
      dictDir,
      U.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      )
     )
  {
    if(! globalParamSet_)
      FatalErrorIn("chPressureGrad::chPressureGrad()")
        << "Must first set global parameters by calling " 
        << "chPressureGrad::initPressureGrad(transportProperties)" << endl
        << abort(FatalError);

    if(chFlowMode_)
    {
        if(specifiedQuantity_ == "gradPbar")
        {
            value_ = mag(gradPbar_);
        }
        else if(specifiedQuantity_ == "varyingGradP")
        {
            value_ = mag(varyingGradP_);
        }
        else if(!gradPDict_.headerOk())
        {
            gradPDict_.add(name, value_);
        }
        else
        {
            gradPDict_.readIfPresent(name, value_);
        }
        
        Info << "Initializing " << name_ << " (from file or by default) as: "
             << name << " = " << value_.value() << endl;
    }
    else
    {
        gradPDict_.writeOpt() = IOobject::NO_WRITE; 
    }

  }
  
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void chPressureGrad::adjust(const volScalarField & rUA)
  {
    if(chFlowMode_)
      {
          volScalarField beta
          (
              IOobject
              (
                  "beta",
                  U_.time().timeName(),
                  U_.mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              scalar(1) - alpha_
          );

          beta.correctBoundaryConditions();

          dimensionedScalar magUbarStar =
              (flowDirection_ & U_)().weightedAverage(beta*alpha_.mesh().V());

          dimensionedScalar volAveAlpha =
              alpha_.weightedAverage(alpha_.mesh().V());

          if(specifiedQuantity_ == "Ubar")
          {
              // Calculate the pressure gradient increment needed to
              // adjust the average flow-rate to the correct value;
              // magUbar: target; magUbarStar: actual in simulation.
              dimensionedScalar gradPplus =
                  (magUbar_ - magUbarStar)
                 /rUA.weightedAverage(U_.mesh().V());

              U_ += flowDirection_*rUA*gradPplus;
              value_ += gradPplus;
              gradPDict_.set(name_, value_);
              
              Info<< solverName_ 
                  << " uncorrected Ubar = " << magUbarStar.value() 
                  << "  "
                  << "with averaged solid volume fraction = "
                  << volAveAlpha.value() << "  "
                  << name_ << " = " << value_.value() << endl;

          }
          else if(specifiedQuantity_ == "gradPbar")
          {

              // dimensionedScalar magUbarStar =
              //     (flowDirection_ & U_)().weightedAverage(beta*alpha_.mesh().V());

              scalar deltaT = U_.mesh().time().value();

              value_ = mag(gradPbar_) + mag(dpdt_)*deltaT;

              Info << solverName_ 
                   << " current Ubar = " << magUbarStar.value() 
                   << "  "
                   << name_ << " = " << value_.value() << endl;
          }
          else if(specifiedQuantity_ == "varyingGradP")
          {

              // dimensionedScalar magUbarStar =
              //     (flowDirection_ & U_)().weightedAverage(beta*alpha_.mesh().V());

              scalar time = U_.mesh().time().value();
              label nStep = U_.mesh().time().timeIndex();
              scalar PI = Foam::constant::mathematical::pi;

              if(varyingType_ == "sinusoidal")
              {
                value_ = mag(varyingGradP_)*Foam::sin(2*PI*time/(period_.value())+0.5*PI);
              }
              else if(varyingType_ == "square")
              {
                int n = round(time/period_.value()+0.5-SMALL);
                value_ = mag(varyingGradP_)*pow(-1,n);
              }
              else if(varyingType_ == "input")
              {
                Info << "time index is: " << nStep << endl;
                if (nStep <= pressureList_.size())
                {
                    value_ = mag(varyingGradP_)*(pressureList_[nStep - 1]/mag(varyingGradP_).value());
                }
                else
                {
                    value_ = 0*mag(varyingGradP_);
                }
              }
              else
              {
                value_ = 0;
              }

              Info << solverName_ 
                   << " current Ubar = " << magUbarStar.value() 
                   << "  "
                   << name_ << " = " << value_.value() << endl;
          }
      }
  }

  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
