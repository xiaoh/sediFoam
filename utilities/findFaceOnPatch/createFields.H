IOdictionary propsDict
(
    IOobject
    (
        "surfaceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

word surfaceName(propsDict.lookup("surfaceName"));
List<vector> startPointList(propsDict.lookup("startPoint"));
List<vector> endPointList(propsDict.lookup("endPoint"));
