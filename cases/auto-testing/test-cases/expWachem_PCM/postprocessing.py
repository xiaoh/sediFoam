#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

cmd = "cp postProcessing/probes/0/p data/p.dat"
os.system(cmd)

cmd = "cp probes/0/p data/p.dat"
os.system(cmd)

time = [];
pressureLocal = [];
pressureDrop = [];
motionData='data/p.dat';
fData = open(motionData, 'r');
for line in fData.readlines()[5:]:
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    pressureDrop.append(float(data[1])-float(data[2]));
    pressureLocal.append(float(data[1]));
    time.append(data[0]);

p1, = plt.plot(time,pressureLocal,'r-',linewidth = 2)
p2, = plt.plot(time,pressureDrop,'ko',linewidth = 2)
lg = plt.legend([p1,p2],["local pressure","relative pressure"],loc=1)
lg.draw_frame(False)
plt.savefig('data/pressureSignal.pdf');
