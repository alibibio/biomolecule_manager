def reverse(seqs):
    reversed_seqs = []
    for seq in seqs:
        reversed_seqs.append(seq[::-1])
    if len(reversed_seqs) == 1:
        return reversed_seqs[0]

    return reversed_seqs


def complement(seqs, complement_nucleotides=None):
    if complement_nucleotides is None:
        complement_nucleotides = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_seqs = []
    if is_RNA(seqs) == 1:
        print('Please, check that you entered DNA sequence')
        return
    else:
        for seq in seqs:
            complement_seq = ''
            for nucleotide in seq:
                complement_seq += complement_nucleotides.get(nucleotide)
            complement_seqs.append(complement_seq)

            if len(complement_seqs) == 1:
                return complement_seqs[0]

            return complement_seqs


def transcribe(seqs):
    if is_RNA(seqs) == 1:
        transcribed_chains = []
        for seq in seqs:
            transcribed_chain = seq.replace('U', 'T').replace('u', 't')
            transcribed_chains.append(transcribed_chain)
        print('Attention! Reversed transcription was done too')

    else:
        transcribed_chains = []
        for seq in seqs:
            transcribed_chain = seq.replace('T', 'U').replace('t', 'u')
            transcribed_chains.append(transcribed_chain)

    if len(transcribed_chains) == 1:
        return transcribed_chains[0]

    return transcribed_chains


def reverse_complement(seqs):
    if is_RNA(seqs) == 1:
        print('Please, check that you entered DNA sequence')
        return
    else:
        complement_list = []
        if type(complement(seqs)) == str:
            complement_list.append(complement(seqs))
        else:
            complement_list = complement(seqs)
        reverse(complement_list)

        return reverse(complement_list)


def is_nucleic_acid(seqs, alphabet=None):
    if alphabet is None:
        alphabet = {'A', 'T', 'G', 'C', 'U'}
    false_amount = 0
    for seq in seqs:
        seq = seq.upper()
        unique_chars = set(seq)
        false_amount += 1 - (unique_chars <= alphabet)

    if false_amount == 0:
        return True
    else:
        raise ValueError('Please, check the correctness of the sequences')


def is_RNA(seqs):
    for seq in seqs:
        seq = seq.upper()
        if 'U' in seq:
            return True


def run_dna_rna_tools(*arguments):
    command = arguments[-1].lower()
    seqs = arguments[:-1]

    if (is_nucleic_acid(seqs)) == 1:
        if command == 'transcribe':
            return transcribe(seqs)
        if command == 'reverse':
            return reverse(seqs)
        if command == 'complement':
            return complement(seqs)
        if command == 'reverse_complement':
            return reverse_complement(seqs)
        else:
            print('Please check the correctness of the command')
