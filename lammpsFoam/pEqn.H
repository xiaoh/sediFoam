volScalarField rUbA = 1.0/UbEqn.A()*beta;

for (int corr = 0; corr < nCorr; corr++)
{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = scalar(1) - alphaf;

    // volScalarField rUaA = 1.0/UaEqn.A();

    // rUbA = 1.0/UbEqn.A()*beta;

    // surfaceScalarField rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf = fvc::interpolate(rUbA);

    // Ua = rUaA*UaEqn.H();
    Ub = rUbA*UbEqn.H()/beta;

    // surfaceScalarField phiDraga =
    //    fvc::interpolate(beta/rhoa*dragCoef*rUaA)*phib
    //  + rUaAf*(g & mesh.Sf());
    surfaceScalarField phiDragb =
        fvc::interpolate(rUbA/rhob)*(fvc::interpolate(cloud.Asrc()) & mesh.Sf())
      + rUbAf*(g & mesh.Sf());

    // fvc::interpolate(dragCoef/rhob*rUbA)*phia + rUbAf*(g & mesh.Sf());


    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            // phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }

    Ua.correctBoundaryConditions();

    phia = (fvc::interpolate(Ua) & mesh.Sf());
      //  + fvc::ddtPhiCorr(rUaA, Ua, phia)
      //  + phiDraga;

    phib =
        (fvc::interpolate(Ub) & mesh.Sf())
      + rUbAf*fvc::ddtCorr(Ub, phib)
      + phiDragb;

    phi = alphaf*phia + betaf*phib;


    surfaceScalarField Dp("(rho*(1|A(U)))", betaf*rUbAf/rhob);

    // Info<< "nNonOrthCorr: " << nNonOrthCorr << endl;
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();
        // Not completely clear about how the nonOrth correction
        // should be modified. Need more thinking.
        if (nonOrth == nNonOrthCorr)
        {
            // Info<< "Dp: " << Dp << endl;
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            // phia -= rUaAf*SfGradp/rhoa;
            phib -= rUbAf*SfGradp/rhob;
            phi = alphaf*phia + betaf*phib;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

            // Ua += (fvc::reconstruct(phiDraga - rUaAf*SfGradp/rhoa));
            // // Ua += rUaA*(fvc::reconstruct(phiDraga/rUaAf - SfGradp/rhoa));
            // Ua.correctBoundaryConditions();

            Ub += (fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob));
            //Ub += rUbA*(fvc::reconstruct(phiDragb/rUbAf - SfGradp/rhob));
            Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

gradP.adjust(rUbA);

Ub.correctBoundaryConditions();

volTensorField gradUb("gradUb",fvc::grad(Ub));
volScalarField nuEff("nuEff",continuousPhaseTurbulence->nuEff());
volScalarField k("k",continuousPhaseTurbulence->k());

B = ((2.0/3.0)*I)*k - nuEff*(twoSymm(gradUb));
