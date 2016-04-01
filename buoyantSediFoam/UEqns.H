// fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);

betaf = fvc::interpolate(beta);
betaPhib = betaf*phib;

{
    UbEqn =
    (
        fvm::ddt(beta,Ub)
      + fvm::div(betaPhib, Ub, "div(phib,Ub)")
      - fvm::Sp(fvc::ddt(beta) + fvc::div(betaPhib), Ub)
      + (Cvm*rhob*alpha*beta/rhob)*
        (
            fvm::ddt(Ub)
          + fvm::div(phib, Ub, "div(phib,Ub)")
          - fvm::Sp(fvc::div(phib), Ub)
        )

        // divDevReff(U) = - laplacian(beta*U)
      + continuousPhaseTurbulence->divDevReff(Ub)
      + continuousPhaseTurbulence->nuEff()*(fvc::grad(beta) & fvc::grad(Ub))
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - beta*fvm::Sp(dragCoef/rhob, Ub)
    //+ alpha/rhob*dragCoef*Ua   // Explicit drag transfered to p-equation
      + beta*alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
      + fvc::average(beta)*gradP.flowDirection()*gradP.value() // Adding pressure gradient
    );

    if (addIBMForce)
    {
        UbEqn -= fvm::Sp(-ibmIndicatorPtr()/ibmRelaxTime, Ub);
    }

    // // Set the velocity directly
    // if (addIBMForce)
    // {
    //     vectorField zerosVelField(ibmIndicatorList.size(), Foam::vector::zero);
    //     UbEqn.setValues(ibmIndicatorList, zerosVelField);
    // }

    UbEqn.relax();

}
