import numpy
import matplotlib
import matplotlib.pyplot as plt
from numpy.fft import *
from pylab import *
import wave
import array
from Bio import SeqIO
import os
import argparse

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

    l = int(len(xhatA))
    xhat=[(xhatA[i]**2+xhatT[i]**2+xhatC[i]**2+xhatG[i]**2)*2/len(seq) for i in range(int(l/2),l)]
    return(xhat,[float(i)/len(xhat) for i in range(0,len(xhat))])


def P(spectrum):
    x=len(spectrum)
    k=int(min(x/3,50))
    peak=max(spectrum[int(x/3)-k:int(x/3)+k])
    p=x*peak/sum(spectrum)
    return(p)


def main(inputfile, outputdir):
    records = list(SeqIO.parse(inputfile, "fasta"))
    for record in records: 
        outfigurename = record.id+'.png'
        outfilename = os.path.join(outputdir, outfigurename)

        plt.gcf()
        ff_p, freq_p=fourier_genome(str(record.seq))
        p=P(ff_p)

        plt.figure(1)
        plt.plot(freq_p,ff_p)
        plt.title("Fourier Transform of %s" % record.id)
        plt.xlabel("frequency")
        plt.ylabel("S")
        plt.text(0.75,15,"P=%.2f"%p)

        plt.savefig(outfilename)
        plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='input fasta file')
    parser.add_argument('-o', '--outputdir', required=True, help='Output directory')
    args = parser.parse_args()

    main(args.input, args.outputdir)