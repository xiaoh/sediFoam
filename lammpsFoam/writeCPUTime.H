Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << "  OpenFOAM/evolve/calcTcField/diffusion/particle move = (" 
    << splitTime[0] << ", "
    << splitTime[1] << ", "
    << splitTime[2] << ", "
    << cloud.diffusionTimeCount()[0] << ", "
    << cloud.diffusionTimeCount()[1] << ", "
    << cloud.particleMoveTime() << ") s"
    << nl << endl;

Info<< "assemble/transpose/flatten/foam->lammps/lammps/lammps->foam = (" 
    << cloud.cpuTimeSplit()[0] << ", "
    << cloud.cpuTimeSplit()[1] << ", "
    << cloud.cpuTimeSplit()[2] << ", "
    << cloud.cpuTimeSplit()[3] << ", "
    << cloud.cpuTimeSplit()[4] << ", "
    << cloud.cpuTimeSplit()[5] << ") s"
    << nl << endl;
