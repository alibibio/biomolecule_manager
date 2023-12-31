import os

def is_suitable_r_quality(seq_quality: str, quality_thresold: int) -> bool:
    """
    Checks if mean quality of sequence is suitable for specified quality
    Arguments:
        seq_quality (str) - string of each nucleotide quality in ASCII encoding
        quality_thresold (int) - the quality, which we admit
    Return:
        True (bool) - mean quality is suitable
        False (bool) - mean quality is not suitable
    """

    qual_summ = 0
    for nucltd_q in seq_quality:
        qual_summ += ord(nucltd_q) - 33
    mean_qual = qual_summ / (len(seq_quality))
    if mean_qual > quality_thresold:
        return True
    else:
        return False


def is_suitable_gc_cont(seqs_seq: str, gc_bound_min: int, gc_bound_max: int) -> bool:
    """
    Checks whether the gc-content of sequence corresponds to a given interval
    Arguments:
        seqs_seq (str) - string of given sequence
        gc_bound_min (int) - the minimum of admited gc-content
        gc_bound_max (inr) - the maximum of admited gc-content
    Return:
        True (bool) - gc-content is suitable
        False (bool) - gc-content is not suitable
    """
    gc_content = (seqs_seq.count('G') + seqs_seq.count('С')) / len(seqs_seq) * 100
    if gc_bound_min <= gc_content < gc_bound_max:
        return True
    else:
        return False


def is_suitable_lenght(seqs_seq: str, length_bound_min: int, length_bound_max: int) -> bool:
    """
    Checks whether the length of sequence corresponds to a given interval
    Arguments:
        seqs_seq (str) - string of given sequence
        length_bound_min (int) - the minimum of admited length
        length_bound_max (inr) - the maximum of admited length
    Return:
        True (bool) - length is suitable
        False (bool) - length is not suitable
    """
    if length_bound_min <= len(seqs_seq) < length_bound_max:
        return True
    else:
        return False


def fastq_tools(seqs, gc_bounds=(0, 101), length_bounds=(0, 2 ** 32), quality_thresold=0) -> set:
    """
    Filters the sequences suitable for given bounds
    Arguments:
        seqs (set) - set of given fastq parameters; {'name of read': ('nucleic acid sequence', 'corresponding quality of read')}
        gc_bounds (int, tuple) - allowed gc-content
        length_bounds (int, tuple) - allowed length
        quality_thresold (int) - allowed quality
    Return:
        filtered_seqs - set of suitable sequences {'name of read': ('nucleic acid sequence', 'corresponding quality of read')}
    """

    if type(gc_bounds) == int:
        gc_bound_max = gc_bounds
        gc_bound_min = 0
    elif type(gc_bounds) == tuple:
        gc_bound_max = gc_bounds[1]
        gc_bound_min = gc_bounds[0]

    if type(length_bounds) == int:
        length_bounds_min = 0
        length_bounds_max = length_bounds
    if type(length_bounds) == tuple:
        length_bounds_min = length_bounds[0]
        length_bounds_max = length_bounds[1]

    filter_counter = dict.fromkeys(seqs.keys(), 0)
    for seqs_nums in seqs.keys():
        seqs_seq = seqs[seqs_nums][0]
        seq_quality = seqs[seqs_nums][1]

        # quality checking
        if is_suitable_r_quality(seq_quality, quality_thresold):
            filter_counter[seqs_nums] += 1
        # gc-content checking
        if is_suitable_gc_cont(seqs_seq, gc_bound_min, gc_bound_max):
            filter_counter[seqs_nums] += 1
        # length checking
        if is_suitable_lenght(seqs_seq, length_bounds_min, length_bounds_max):
            filter_counter[seqs_nums] += 1

    filtered_seqs = {}
    for f_count in filter_counter.keys():
        if filter_counter.get(f_count) == 3:
            filtered_seqs[f_count] = seqs[f_count]

    return filtered_seqs


def read_fastq(input_filename : str) -> dict:
    """
    Makes dictionary of fastq strings values
    Args:
        input_filename (str): the name of fastq file with extension .fastq
    Returns:
        dictionary of fastq values in format: {name: (sequence, quality)}
    """
    counter = 3
    samples = []
    with open(input_filename) as fastq_file:
        sample = []

        for line in fastq_file:
            if not line.startswith('+'):
                sample.append(line.strip())
            if len(sample) == counter:
                samples.append(sample)
                sample = []

    fastq_dict = {}
    for sample in samples:
        fastq_dict[sample[0].split()[0]] = tuple(sample[1:])

    return fastq_dict


def write_fastq(fastq_output_dict: dict, output_filename='example_output.fastq', output_folder='fastq_filtrator_resuls'):
    """
    Args:
        fastq_output_dict(dict): dictionary of fastq values in format: {name: (sequence, quality)}
        output_filename(str): name of output file. Default value 'example_output.fastq'
        output_folder(str): name of folder with results. Default value 'fastq_filtrator_resuls'

    Returns:
        Creates if it necessary new folder with file .fastq with given reads
    """

    if not os.path.exists(output_folder):
        os.makedirs('fastq_filtrator_resuls')

    with open(os.path.join(output_folder, output_filename), mode='w') as fastq_file:
        for sample in fastq_output_dict:
            fastq_file.write(sample + '\n')
            fastq_file.write(fastq_output_dict[sample][0] + '\n')
            fastq_file.write('+' + sample + '\n')
            fastq_file.write(fastq_output_dict[sample][1] + '\n')
