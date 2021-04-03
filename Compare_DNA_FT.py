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

from DNA_FT import *



def main(inputfile1, inputfile2, outputdir):
    record1 = list(SeqIO.parse(inputfile1, "fasta"))[0]
    record2 = list(SeqIO.parse(inputfile2, "fasta"))[0]

    # compute Fourier transform for each sequence
    ff_p1, freq_p1=fourier_genome(str(record1.seq))
    ff_p2, freq_p2=fourier_genome(str(record2.seq))

    #plot on the same axes
    outfilename = os.path.join(outputdir, 'Comparison_FT.png')
    plt.figure(1)
    plt.plot(freq_p1,ff_p1, label = record1.id)
    plt.plot(freq_p2,ff_p2, label = record2.id)
    plt.legend()
    plt.title("Fourier Transform")
    plt.xlabel("frequency")
    plt.ylabel("S")
    plt.savefig(outfilename)
    plt.clf()

    common_x = set(freq_p1).intersection(set(freq_p2))
    print(len(freq_p1), len(freq_p2))
    #compare spectra using euclidean distance



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', '--input1', required=True, help='input fasta file 1')
    parser.add_argument('-f2', '--input2', required=True, help='input fasta file 2')
    parser.add_argument('-o', '--outputdir', required=True, help='Output directory')
    args = parser.parse_args()

    main(args.input1, args.input2, args.outputdir)