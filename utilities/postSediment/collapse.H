    fileName path(alpha.rootPath()/alpha.caseName()/"graphs"/alpha.instance());
    mkDir(path);

    scalarField pValues(channelIndexing.collapse(p));
    scalarField alphaValues(channelIndexing.collapse(alpha));

    const scalarField& y = channelIndexing.y();

    makeGraph(y, pValues, "p", path, gFormat);
    makeGraph(y, alphaValues, "alpha", path, gFormat);
