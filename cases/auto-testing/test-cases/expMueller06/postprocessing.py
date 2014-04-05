#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt

cmd = "sample -latestTime"
os.system(cmd)

cmd = "cp postProcessing/sets/*/* data/"
os.system(cmd)

cmd = "cp sets/*/* data/"
os.system(cmd)

for i in [1]:
    plt.figure(i)
    x = [];
    u_y = [];
    x_bench = [];
    u_y_bench = [];
    x_exp = [];
    u_y_exp = [];
    p = [];
    
    pData = 'data/lineY' + str(i+2) + '_UaMean.xy';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x.append(data[0]);
        u_y.append(data[2]);

    p.append(plt.plot(x, u_y, 'r-', linewidth = 4))

    pData = 'data/sets_bench/lineY' + str(i+2) + '_UaMean.xy';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x_bench.append(data[0]);
        u_y_bench.append(data[2]);

    p.append(plt.plot(x_bench, u_y_bench, 'b-.', linewidth = 4))

    pData = 'data/experimentData/Uy' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x_exp.append(0.022+float(data[0]));
        u_y_exp.append(data[1]);

    p.append(plt.plot(x_exp, u_y_exp, 'ko', markersize = 8))

    lg = plt.legend([p[2],p[1],p[0]],\
            ["experiment","previous simulation","present simulation"],loc=4)
    plt.ylim(-0.2*i,0.2*i)
    lg.draw_frame(False)
    pdfName = 'data/timeAveragedVelocity06_' + str(i) + '.pdf'
    plt.savefig(pdfName);
    


for i in [1,2]:
    plt.figure(i+3)
    x = [];
    alpha = [];
    x_bench = [];
    alpha_bench = [];
    x_exp = [];
    alpha_exp = [];
    p = [];
    
    pData = 'data/lineY' + str(i) + '_alphaMean.xy';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x.append(data[0]);
        alpha.append(1-float(data[1]));

    p.append(plt.plot(x, alpha, 'r-', linewidth = 4))

    pData = 'data/sets_bench/lineY' + str(i) + '_alphaMean.xy';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x_bench.append(data[0]);
        alpha_bench.append(1-float(data[1]));

    p.append(plt.plot(x_bench, alpha_bench, 'b-.', linewidth = 4))

    pData = 'data/experimentData/concentration' + str(i) + '.dat';
    fData = open(pData, 'r');
    for line in fData.readlines():
        data = [xx.strip() for xx in line.split(None)]
        if (not data):
            continue
    
        x_exp.append(data[0]);
        alpha_exp.append(data[1]);

    p.append(plt.plot(x_exp, alpha_exp, 'ko', markersize = 8))

    lg = plt.legend([p[2],p[1],p[0]],\
            ["experiment","previous simulation","present simulation"],loc=4)
    plt.ylim(0.2,0.8)
    lg.draw_frame(False)
    pdfName = 'data/fluidVolumeFraction06_' + str(i) + '.pdf'
    plt.savefig(pdfName);
