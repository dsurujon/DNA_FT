import numpy
import matplotlib
import matplotlib.pyplot as plt
from numpy.fft import *
from pylab import *
import wave
import array

filename="PEX3 genomic.txt"
f=open(filename) 
lines=f.readlines()
f.close()
pex3="".join(lines[1:])
pex3=pex3.replace("\n","")
#if the file has intergenic regions 1000bp upstream and 1000bp downstream 
pex3=pex3[1000:-1000]
def fourier_genome(seq):
    xA=[0 for i in range(0,len(seq))]
    xT=[0 for i in range(0,len(seq))]
    xG=[0 for i in range(0,len(seq))]
    xC=[0 for i in range(0,len(seq))]

    for i in range(0,len(seq)):
        if seq[i]=="A":xA[i]=1
        if seq[i]=="T":xT[i]=1
        if seq[i]=="G":xG[i]=1
        if seq[i]=="C":xC[i]=1

    xhatA=abs(fft(xA))
    xhatG=abs(fft(xG))
    xhatT=abs(fft(xT))
    xhatC=abs(fft(xC))
    xhat=[(xhatA[i]**2+xhatT[i]**2+xhatC[i]**2+xhatG[i]**2)*2/len(seq) for i in range(len(xhatA)/2,len(xhatA))]
    #xhat=[(xhatA[i])*2/len(seq) for i in range(len(xhatA)/2,len(xhatA))]
    return xhat,[float(i)/len(xhat) for i in range(0,len(xhat))]
def P(spectrum):
    x=len(spectrum)
    k=min(x/3,50)
    peak=max(spectrum[(x/3)-k:(x/3)+k])
    p=x*peak/sum(spectrum)
    return p

ff_p, freq_p=fourier_genome(pex3)
p=P(ff_p)

plt.figure(1)
plt.plot(freq_p,ff_p)
plt.title("Fourier Transform of the pex3 locus")
plt.xlabel("frequency")
plt.ylabel("S")
plt.text(0.75,15,"P=%.2f"%p)

show()
