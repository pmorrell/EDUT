# EDUT Basic Usage

 example command line execution: `./EDUT.pl -i a_fasta_file.txt -b -d none -s`

1. Make the EDUT.pl file executable with `chmod +x EDUT.pl`.

2. Create an aligned fasta file. The file should be formatted so that each sequence name and each sequence is a single line of text with no line breaks within a sequence. Line breaks follow each sequence name and each nucleotide sequence, including the last sequence. Use UNIX line breaks. All sequences in the data set must be the same length.

3. The program can be executed with command line options. The first option which is required for evaluating DNA sequence data is -i followed by the name of the fasta file.

4. You can also use simulated data from Hudson's simulation tool mksample (ms) as input. This allows you to see the expected results of the error detection in the absence of phasing or sequencing errors. Simulations can be generated with your estimated population parameters. In this way for example you can see how much homologous gene conversion would elevate the weighted priority scores for your data. To use this feature you should use the option -a followed by the name of the simulation file. An example of the format for the simulated data file is provided. The -a option is to be used instead of -i. When using simulated data as input the options that follow can still be used, including changing the delimiter for the simulation binary data.

5.  The program can output your data in an ms style binary format. This format is used by other programs and provides a summary of the haplotypes. Often this binary format allows the user to visually assess the influence of recombination.
		
	-b indicates print binary data
	
	-d followed by one of the following words (typed without the quotes) for the desired delimiter for the binary data output:
	
	"tab" is the default delimiter.
	"space"
	"none"

6. Choose the type of output for the error detection function.
	
	-v indicates verbose output. This output lists all of the triplets in the described pattern and reports the counts for all four configurations for each of the two pairs (AB and BC) making up the triplet.
	
    -m indicates a moderate amount of output. This will output the pairs that have been flagged as unlikely configurations.
    
	-s indicates short summary output. This is the recommended option because it reports the raw and a weighted priority score. The summary output can also be combined with the verbose or moderate output.
	

