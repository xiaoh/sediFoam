{
    word scheme("div(phi,alpha)");

    surfaceScalarField phir = phia - phib;

    Info<< "Max Ur Courant Number = "
        << (
               max
               (
                   mesh.surfaceInterpolation::deltaCoeffs()*mag(phir)
                  /mesh.magSf()
                )*runTime.deltaT()
            ).value()
        << endl;

    // for (int acorr = 0; acorr < nAlphaCorr; acorr++)
    {
      /*
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phi, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, beta, scheme), alpha, scheme)
        );
        alphaEqn.relax();
        alphaEqn.solve();
      */

      // Set volume fraction from Lammps computation.
      // A fixed value set in initial condition for now.
      // No need to change here.
      // alpha =

        /*
        fvScalarMatrix betaEqn
        (
              fvm::ddt(beta)
            + fvm::div(phi, beta, scheme)
            + fvm::div(-fvc::flux(phir, scalar(1) - beta, scheme), beta, scheme)
        );
        betaEqn.relax();
        betaEqn.solve();

        alpha =
            0.5*(scalar(1) + sqr(scalar(1) - beta)
          - sqr(scalar(1) - alpha));
        */
        beta = scalar(1) - alpha;
        alpha.correctBoundaryConditions();
        beta.correctBoundaryConditions();
    }

    Info<< "Dispersed phase volume fraction = "
        << alpha.weightedAverage(mesh.V()).value()
        << "  Min(alpha) = " << min(alpha).value()
        << "  Max(alpha) = " << max(alpha).value()
        << endl;
}
