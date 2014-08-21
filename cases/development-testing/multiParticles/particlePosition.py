#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

p = [];
for i in [1,2,3]:
    cmd = 'grep "^' + str(i) +' 0.00" snapshot.bubblemd > data/p' + str(i) + '.dat';
    os.system(cmd);
    x_p = [];
    y_p = [];
    x_p_bench = [];
    y_p_bench = [];
    
    pData='data/p' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [x.strip() for x in line.split(None)]
        if (not data):
            continue
    
        x_p.append(data[4]);
        y_p.append(data[5]);

    p.append(plt.plot(x_p,y_p,'k-o',markersize = 5))

    pData='data/origin/p' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [x.strip() for x in line.split(None)]
        if (not data):
            continue
    
        x_p_bench.append(data[4]);
        y_p_bench.append(data[5]);

    p.append(plt.plot(x_p_bench,y_p_bench,'r-o',markersize = 3))

lg = plt.legend([p[0],p[1]],["current result","benchmark"],loc=4)
lg.draw_frame(False)
plt.savefig('data/multiParticlesPositionDia.pdf');
