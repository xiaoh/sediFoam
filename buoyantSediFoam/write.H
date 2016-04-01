// if (runTime.outputTime())
{
    volVectorField Ur
    (
        IOobject
        (
            "Ur",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua - Ub
    );

    runTime.write();
}
