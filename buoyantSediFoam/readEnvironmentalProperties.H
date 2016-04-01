Info<< "\nReading environmentalProperties" << endl;

IOdictionary environmentalProperties
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedVector g(environmentalProperties.lookup("g"));
