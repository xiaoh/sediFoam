{
    volVectorField Ur = Ua - Ub;
    volScalarField magUr = mag(Ur);

    // scalar maxPossibleAlpha = scalar(0.70);
    if (max(alpha).value() > maxPossibleAlpha)
    {
        Info<< "Warning: Unphysical large solid volume fraction "
            << max(alpha).value() << endl;
        Info<< "A cap of " << maxPossibleAlpha << " is applied." << endl;

        alpha = Foam::min(alpha, maxPossibleAlpha);
        beta = scalar(1) - alpha;
    }

cloud.calcTcFields();

dragCoef = cloud.Omega();

// Lift force not considered for now.
// However, rhoa is eliminated from the formula
// (why should the lift force depends on the solid density?)
liftCoeff = Cl*beta*rhob*(Ur ^ fvc::curl(U));
// liftCoeff = Cl*(beta*rhob + alpha*rhoa)*(Ur ^ fvc::curl(U));
}
