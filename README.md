## Fourier spectra of DNA sequences
The `DNA_FT.py` script will read all sequences in the input fasta file (in this case `test.fasta`), and generate the Fourier spectrum plot for each fasta entry. The plots for each entry are saved as `.png` files in the output directory specified with `-o`. 

All dependencies are included in the `dna_ft.yml` file.     
Build and activate the conda environment: 
```
conda env create -f dna_ft.yml 
conda activate dna_ft
```
Then run the `DNA_FT.py` script:
```
python python DNA_FT.py -i test.fasta -o ./output/
```

### Comparing two sequences
It's possible to generate an overlapping plot for two fasta files (each with a single record - if there are multiple records, only the first will be considered from each file). For this, use the `Compare_DNA_FT.py` script: 
```
python Compare_DNA_FT.py -f1 NC_045512.2.fasta -f2 MT066176.1.fasta -o output
```

The code uses the method outlined in the following publication:    

Tiwari, S., Ramachandran, S., Bhattacharya, A., Bhattacharya, S. & Ramaswamy, R. Prediction of probable genes by Fourier analysis of genomic sequences. Comput. Appl. Biosci. 13, 263–270 (1997).
