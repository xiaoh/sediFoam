#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

linestyle = ['k-','ko','b-','bo','r-','ro','g-','go']
p = [];
for i in [1,2,3]:
    cmd = 'grep "^' + str(i) +' 0.00" snapshot.bubblemd > data/p' + str(i) + '.dat';
    os.system(cmd);
    t_now = 0;
    timeStep = 0.0002;
    time = [];
    u_p = [];
    u_bench = [];
    
    pData='data/p' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [x.strip() for x in line.split(None)]
        if (not data):
            continue
    
        u_p.append(data[8]);
        time.append(t_now);
        t_now = t_now+timeStep;

    p.append(plt.plot(time,u_p,linestyle[2*i-2],linewidth = 3))

    pData='data/origin/p' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [x.strip() for x in line.split(None)]
        if (not data):
            continue
    
        u_bench.append(data[8]);

    p.append(plt.plot(time,u_bench,linestyle[2*i-1],markersize = 6))

lg = plt.legend([p[0],p[1],p[2],p[3],p[4],p[5]],\
        ["particle #1 (current)","particle #1 (benchmark)",\
         "particle #2 (current)","particle #2 (benchmark)",\
         "particle #3 (current)","particle #3 (benchmark)"],\
          loc=1)
lg.draw_frame(False)
plt.savefig('data/multiParticlesVelocityDia.pdf');
