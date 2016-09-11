    // adding the force based on immerse boundary method 
    // can be used when a constant flow velocity is need in some region
    Switch addIBMForce(transportProperties.lookupOrDefault<Switch>("addIBMForce", false));

    if (addIBMForce)
    {
        Info<< "NOTICE!! Adding IBM force in present simulations." << endl;
    }

    autoPtr<volScalarField> ibmIndicatorPtr;

    dimensionedScalar ibmRelaxTime
    (
        dimensioned<scalar>::lookupOrDefault
        (
            "ibmRelaxTime",
            transportProperties,
            runTime.deltaT().value()*scalar(3),
            dimTime
        )
    );

    DynamicList<label> ibmIndicatorList(0);

    if (addIBMForce)
    {
        ibmIndicatorPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "ibmIndicator",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        volScalarField& ibmIndicatorField = ibmIndicatorPtr();

        forAll(ibmIndicatorField, cellI)
        {
            if (ibmIndicatorField[cellI] > 0)
            {
                ibmIndicatorList.append(cellI);
            }
        }
    }
