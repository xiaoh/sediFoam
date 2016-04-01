Info<< "Moving Particles" << endl;

cloud.evolve();
cloud.averageInfo();

// Update Uc field as well
Uc = Ua*alpha;
