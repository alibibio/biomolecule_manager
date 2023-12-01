import src.dna_rna_tools as dna_rna
import src.fastq_tools as fastq
import src.protein_tools as protein


def dna_rna_tools(*args):
    return dna_rna.run_dna_rna_tools(*args)


def fastq_tools(input_path: str, gc_bounds=(0, 101), length_bounds=(0, 2 ** 32), quality_thresold=0):
    seqs = fastq.read_fastq(input_path)
    fastq_output = fastq.fastq_tools(seqs, gc_bounds=gc_bounds, length_bounds=length_bounds, quality_thresold=quality_thresold)
    fastq.write_fastq(fastq_output)


def protein_tools(function: str, *prots: str):
    return protein.protein_tools(function, *prots)
