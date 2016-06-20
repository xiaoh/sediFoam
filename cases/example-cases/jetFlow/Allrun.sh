blockMesh > log.blockMesh
decomposePar > log.decomposePar

mpirun -np 6 lammpsFoam -parallel > log.parallel
