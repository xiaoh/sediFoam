Info<< "Moving Particles" << endl;

cloud.evolve();

// Update Uc field as well
Uc = Ua*alpha;
