#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

cmd = "cp postProcessing/probes/0/p data/p.dat"
os.system(cmd)

cmd = "cp probes/0/p data/p.dat"
os.system(cmd)

time = [];
pressureDrop = [];
motionData='data/p.dat';
fData = open(motionData, 'r');
for line in fData.readlines()[4:]:
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    pressureDrop.append(float(data[1])-float(data[2]));
    time.append(data[0]);

time_bench = [];
pressureDrop_bench = [];
motionData='data/p_bench.dat';
fData = open(motionData, 'r');
for line in fData.readlines():
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    pressureDrop_bench.append(data[1]);
    time_bench.append(data[0]);

p1, = plt.plot(time,pressureDrop,'r-',linewidth = 4)
p2, = plt.plot(time_bench,pressureDrop_bench,'ko',markersize = 8)
lg = plt.legend([p2,p1],["benchmark","current code"],loc=1)
lg.draw_frame(False)
plt.savefig('data/pressureDrop.pdf');
