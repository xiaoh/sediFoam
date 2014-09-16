// fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);

{
    UbEqn =
    (
        (scalar(1) + Cvm*rhob*alpha/rhob)*
        (
            fvm::ddt(Ub)
          + fvm::div(phib, Ub, "div(phib,Ub)")
          - fvm::Sp(fvc::div(phib), Ub)
        )

        // divDevReff(U) = - laplacian(beta*U)
      + 1/beta*continuousPhaseTurbulence->divDevReff(Ub) 
      + 1/beta*continuousPhaseTurbulence->nuEff()*(fvc::grad(beta) & fvc::grad(Ub))
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - fvm::Sp(dragCoef/rhob, Ub)
    //+ alpha/rhob*dragCoef*Ua   // Explicit drag transfered to p-equation
      + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
      + gradP.flowDirection()*gradP.value() // Adding pressure gradient
    );

    UbEqn.relax();
}
