    IOobject phibHeader
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );


    surfaceScalarField* phibPtr;

    if (phibHeader.headerOk())
    {
        Info<< "Reading face flux field phib\n" << endl;

        phibPtr = new surfaceScalarField
        (
            IOobject
            (
                "phib",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Calculating face flux field phib\n" << endl;

        wordList phiTypes
        (
            Ub.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        for (label i = 0; i < Ub.boundaryField().size(); i++)
        {
            if (isType<fixedValueFvPatchVectorField>(Ub.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phibPtr = new surfaceScalarField
        (
            IOobject
            (
                "phib",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(Ub) & mesh.Sf(),
            phiTypes
        );
    }

    surfaceScalarField& phib = *phibPtr;
