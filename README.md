# motif discovery among the TRN promoters

This contains a collection of codes that pull out the promoter region of TRN genes in C.elegans, and tries to find the common cis-regulatory motif.
In theory, this should work for any list of genes for a given cell fate.

The "TRN promoters.Rmd" is a set of R codes that isolate the promoter sequences (1000 bp or 500 bp upstream) of TRN specific genes (in total 187 genes)

The sequences are stored in "TRN_promoter_sequences.txt" (1000 bp promoter region) with the gene information indexed in "TRN_promoters_granges.csv" (the csv file can be read in as grange files using the Bioconductor packages in R)

500 bp promoter sequences are stored in "TRN_promoter_sequences_500.txt" with the gene information indexed in "TRN_promoters_granges_500.csv"

"GibbsMotifSearch.py" is a set of python codes I wrote to discover the common motif among the TRN promoters. Three numbers ("K","N","r") need to be added to the sequence file as the first line, separated by space. The idea is to randomly choose a Kmer from each promoter sequence to construct a matrix. Score the matrix. Then, randomly choose a Kmer and replace it with another Kmer randomly drawn from the same promoter (not really random, weighed random). Score the matrix again, if the score improves, update the matrix. Repeat this process N times. 

Reconstruct the set of random motifs, redo the above for "r" times. At the end, compare all matrixes and identify the best. 

I did two test runs.
"run_1.txt" contains the input file: K = 15, N = 200, r = 2000, 1000 bp promoters
"run_1_results.txt" shows the results (a set of motifs and the consensus)
"run_2.txt" contains the input file: K = 15, N = 2000, r = 2000, 500 bp promoters
"run_2_results.txt" shows the results (a set of motifs and the consensus)

