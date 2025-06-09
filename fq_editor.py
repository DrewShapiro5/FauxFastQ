import random

class FQEditor:

    def __init__(self, seed: int, file_forward, file_reverse, unidirectional: bool = False, file_mode: str = 'w', preserve_case: bool = False, use_random_bases: bool = False):
        self._header_str = "@M05774:236:000000000-LKC5Y:1:1107:27241:{index} {dir}:N:0:GATCAGAT+TCCGCGAA"
        self._random = random.Random()
        self._random.seed(seed)
        self._unidirectional = unidirectional
        self._file_mode = file_mode
        self._total_reads = 0
        self.file_forward = file_forward
        self.file_reverse = file_reverse
        self._preserve_case = preserve_case
        self._use_random_bases = use_random_bases

    # Creates fastq files with deletions for each base pair
    def create_simulated_fastqs_deletion(self, sequence, quality_forward, quality_reverse, edit_len):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            edit_func=lambda seq, idx, l: seq[:idx] + seq[idx + l:],
                            iter_len=len(sequence) - edit_len + 1,
                            edit_len=edit_len)

    # Creates fastq files with replacements for each base pair
    def create_simulated_fastqs_replacement(self, sequence, quality_forward, quality_reverse, edit_len):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            edit_func=self.pair_swap,
                            iter_len=len(sequence) - edit_len + 1,
                            edit_len=edit_len)

    # Creates fastq files with insertions between each pair of base pairs
    def create_simulated_fastqs_insertion(self, sequence, quality_forward, quality_reverse, edit_len):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            edit_func=self.pair_insertion,
                            iter_len=len(sequence) + 1,
                            edit_len=edit_len)

    def _process_edits(self, sequence, quality_forward, quality_reverse,
                       edit_func, iter_len, edit_len=1):
        """
        Driver method for creating edited FASTQ files.
        Args:
            edit_func: Function that performs the sequence edit (deletion, swap, or insertion)
            edit_len: Length of edit (default=1)
        """
        if not self._preserve_case:
            sequence = sequence.lower()
        # Determine output files based on unidirectional setting
        output_files = []
        if self._unidirectional:
            output_files.append(open(self.file_forward, self._file_mode))
        else:
            output_files.append(open(self.file_forward, self._file_mode))
            output_files.append(open(self.file_reverse, self._file_mode))

        for index in range(iter_len):
            # Generate headers
            header_forward = self._header_str.format(index=index + self._total_reads, dir=1)
            header_reverse = self._header_str.format(index=index + self._total_reads, dir=2)

            # Apply the edit
            edited_sequence = edit_func(sequence, index, edit_len)

            # Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
            sequence_forward = edited_sequence[:len(quality_forward)]
            sequence_reverse = edited_sequence[-len(quality_reverse):][::-1]

            # Write entries
            write_fastq_entry(output_files[0], header_forward,
                                    sequence_forward, quality_forward)

            if not self._unidirectional:
                write_fastq_entry(output_files[1], header_reverse,
                                        dna_complement(sequence_reverse),
                                        quality_reverse)
        # Ensure files are closed
        for f in output_files:
            print("{} edits written to file {}".format(iter_len, f.name))
            f.close()
        self._total_reads += iter_len

    # Replaces the base pair at index 'index' with a random base pair
    def pair_swap(self, sequence, index, edit_len):
        replacement = []
        list_sequence = list(sequence)
        if self._use_random_bases:
            for i in range(edit_len):
                replacement.append(random_base(self._random))
        else:
            replacement = dna_complement(circular_slice(list_sequence, index, edit_len))
        list_sequence = list_sequence[:index] + [replacement.upper()] + list_sequence[index + edit_len:]
        return ''.join(list_sequence)

    # Inserts the complement of the pair at index 'index' after the pair
    def pair_insertion(self, sequence, index, edit_len):
        insertion = ''
        list_sequence = list(sequence)
        if self._use_random_bases:
            for i in range(edit_len):
                insertion += random_base(self._random)
        else:
            insertion = dna_complement(circular_slice(list_sequence, index, edit_len))
        list_sequence.insert(index, ''.join(insertion).upper())
        return ''.join(list_sequence)

def circular_slice(arr, start, length):
    result = []
    for i in range(length):
        index = (start + i) % len(arr)
        result.append(arr[index])
    return result

# Returns the complement of a given sequence
def dna_complement(sequence):
    complement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c'
    }
    complement = []
    for base in sequence:
        complement.append(complement_dict[base])
    return ''.join(complement)

def random_base(r: random.Random):
    base_arr = ['A', 'C', 'G', 'T']
    return base_arr[r.randint(0, len(base_arr) - 1)]

def write_fastq_entry(file, header, sequence, quality):
    """Helper method to write a single FASTQ entry."""
    file.write(f'{header}\n{sequence}\n+\n{quality}\n')
