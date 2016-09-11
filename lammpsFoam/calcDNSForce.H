
if (addDNSForce)
{
    turbulenceForce.internalField() = ReImSum
    (
        fft::reverseTransform
        (
            K/(mag(K) + 1.0e-6) ^ forceGen.newField(), K.nn()
        )
    );

    Info<< "k("
         << runTime.timeName()
         << ") = "
         << 3.0/2.0*average(magSqr(Ub)).value() << endl;

    Info<< "epsilon("
         << runTime.timeName()
         << ") = "
         << (
                 0.5*nub*average
                 (
                     magSqr(fvc::grad(Ub) + fvc::grad(Ub)().T())
                 )
            ).value() << endl;

    Info<< "U.f("
         << runTime.timeName()
         << ") = "
         << 181.0*average(Ub & turbulenceForce).value() << endl;
}
