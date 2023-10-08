# Biomolecule manager
Contains tools for working with proteins, nucleic acids and fastq data sorting

Import biomolecules_manager library, then you can call functions with certain arguments:
- dna_rna_tools
- fastq_tools
- protein_tools\
Read about each function below

## DNA-RNA tools
Provide the function arguments:
1) the first argument should be command (choose one of: transcribe, reverse, complement, reverse_complement)
2) other arguments should be the sequencee(s) of the nucleic acid(s) 

Here is the catalogue of actions the user can choose:

- *transcribe* - makes transcribed RNA or DNA(reversed transcription) chain
- *complement* - returns complement chain(s) for DNA sequences
- *reverse* - makes reversed chain(s)
- *reverse_complement* - gives reversed complement sequense(s) of DNA chain

## Protein tools
Provide the tool with the sequence(s) of the protein(s) in 1-letter format (for example, DYKDDDDK) and the function needed. If you
occasionally write down a non-peptide sequence, the programm will return an error.  

Here is the catalogue of actions the user can choose: 

- *count_length*: gives the length(s) of the protein sequence(s)  
- *count_nucleotide_length*: counts the length(s) of the coding nucleotide sequence(s) of the protein sequence(s)  
- *count_molecular_mass*: calculates molecular mass of the input (the algorithm takes into consideration water mass and subtracts it)    
- *show_content*: shows the aminoacid content of the protein(s)  
- *convert_1_to_3*: converts 1-letter format into 3-letter one  
- *count_extinction_280nm*: counts the molar extinction coefficient (this function counts cystine contribution to extinction coefficient as two cysteins give 1 SS-bond) 


## Fastq tools
 Filters the sequences suitable for given bounds. 
Arguments:
 Accept sequences from fastq files as set. Example:{'name of read': ('nucleic acid sequence', 'corresponding quality of read')}
 Parameters can be set using nominal variables
 
Here is the catalogue of actions the user can set:
*gc_bounds = (y, x)* or *gc_bounds = x*
- filters values between 0 and *x* (doesn't include x) if input is int
- filters values between *y* and *x* (doesn't include x) if input is tuple
- default - (0, 100)

*length_bounds = (y, x)* or *length_bounds = x*
- filters values between 0 and *x* (doesn't include x) if input is int
- filters values between *y* and *x* (doesn't include x) if input is tuple
- default - (0, 2 ** 32)

*quality_thresold* 
- filters values higher then given (int)
- default - 0

Fastq tools returns the filtered set