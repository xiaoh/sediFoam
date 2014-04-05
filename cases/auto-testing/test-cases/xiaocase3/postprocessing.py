#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

cmd = 'grep "1 1" snapshot.bubblemd > data/singleParticle.dat';
os.system(cmd);

varIndex = [0,1,2,3,4,5,6,7,8,9];
dataFloat = [0,0,0,0,0,0,0,0,0,0];
time = [];
velocity = [];

motionData='data/singleParticle.dat';
fData = open(motionData, 'r');
t_now = 0;
timeStep = 0.0002;
for line in fData.readlines():
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    for i in varIndex:
        dataFloat[i] = data[varIndex[i]];

    velocity.append(dataFloat[8]);
    time.append(t_now);
    t_now = t_now+timeStep;

xiao_time = [];
xiao_velocity = [];
xiaoData='data/xiaoCase3.dat';
fData = open(xiaoData, 'r');
for line in fData.readlines():
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    xiao_time.append(data[0]);
    xiao_velocity.append(data[1]);

lammps08_time= [];
lammps08_velocity = [];
lammps08Data='data/lammps08.dat';
fData = open(lammps08Data, 'r');
for line in fData.readlines():
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    lammps08_time.append(data[0]);
    lammps08_velocity.append(data[2]);

p1, = plt.plot(time,velocity,'b-.',linewidth = 4)
p2, = plt.plot(xiao_time,xiao_velocity,'r--',linewidth = 4)
p3, = plt.plot(lammps08_time,lammps08_velocity,'ko',markersize = 10);
lg = plt.legend([p2,p3,p1],["benchmark","lammps08 code","present code"],loc=4)
lg.draw_frame(False)
plt.savefig('data/singleParticle.pdf');
