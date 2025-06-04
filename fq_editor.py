import random

class Fqeditor:

    def __init__(self, seed: int, unidirectional: bool = False):
        self.header_str = "@M05774:236:000000000-LKC5Y:1:1107:27241:{index} {dir}:N:0:GATCAGAT+TCCGCGAA"
        self._random = random.Random()
        self._random.seed(seed)
        self.unidirectional = unidirectional

    # Creates fastq files with deletions for each base pair
    def create_simulated_fastqs_deletion(self, sequence, quality_forward, quality_reverse, edit_len, file_forward, file_reverse):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            file_forward,
                            file_reverse,
                            edit_func=lambda seq, idx, l, r: seq[:idx] + seq[idx + l:],
                            iter_len=len(sequence) - edit_len,
                            edit_len=edit_len)

        """
        with open('single_deletion_forward.fastq', "w") as f1:
            with open('single_deletion_reverse.fastq', "w") as f2:
                for index in range(len(sequence) - edit_len):
                    header_forward = self.header_str.format(index=index, dir=1)
                    header_reverse = self.header_str.format(index=index, dir=2)
                    edited_sequence = sequence[:index] + sequence[index + edit_len:]
                    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                    sequence_forward = edited_sequence[:len(quality_forward)]
                    sequence_reverse = edited_sequence[-len(quality_reverse):][::-1]

                    #Write the contents to the appropriate files
                    f1.write(f'{header_forward}\n{sequence_forward}\n+\n{quality_forward}\n')
                    f2.write(f'{header_reverse}\n{dna_complement(sequence_reverse)}\n+\n{quality_reverse}\n')
            """

    # Creates fastq files with replacements for each base pair
    def create_simulated_fastqs_replacement(self, sequence, quality_forward, quality_reverse, edit_len, file_forward, file_reverse):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            file_forward,
                            file_reverse,
                            edit_func=single_pair_swap,
                            iter_len=len(sequence) - edit_len,
                            edit_len=edit_len)

        """
        with open('single_swap_forward.fastq', "w") as f1:
            with open('single_swap_reverse.fastq', "w") as f2:
                for index in range(len(sequence) - 1):
                    header_forward = self.header_str.format(index=index, dir=1)
                    header_reverse = self.header_str.format(index=index, dir=2)
                    edited_sequence = single_pair_swap(sequence, index)
                    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                    sequence_forward = edited_sequence[:len(quality_forward)]
                    sequence_reverse = edited_sequence[-len(quality_reverse):][::-1]

                    #Write the contents to the appropriate files
                    f1.write(f'{header_forward}\n{sequence_forward}\n+\n{quality_forward}\n')
                    f2.write(f'{header_reverse}\n{dna_complement(sequence_reverse)}\n+\n{quality_reverse}\n')
                    """

    # Creates fastq files with insertions between each pair of base pairs
    def create_simulated_fastqs_insertion(self, sequence, quality_forward, quality_reverse, edit_len, file_forward, file_reverse):
        self._process_edits(sequence,
                            quality_forward,
                            quality_reverse,
                            file_forward,
                            file_reverse,
                            edit_func=single_pair_insertion,
                            iter_len=len(sequence) - edit_len,
                            edit_len=edit_len)

        """
        with open('single_insertion_forward.fastq', "w") as f1:
            with open('single_insertion_reverse.fastq', "w") as f2:
                for index in range(len(sequence) - 1):
                    header_forward = self.header_str.format(index=index, dir=1)
                    header_reverse = self.header_str.format(index=index, dir=2)
                    edited_sequence = single_pair_insertion(sequence, index)
                    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                    sequence_forward = edited_sequence[:len(quality_forward)]
                    sequence_reverse = edited_sequence[-len(quality_reverse):][::-1]

                    #Write the contents to the appropriate files
                    f1.write(f'{header_forward}\n{sequence_forward}\n+\n{quality_forward}\n')
                    f2.write(f'{header_reverse}\n{dna_complement(sequence_reverse)}\n+\n{quality_reverse}\n')
                    """


    def _process_edits(self, sequence, quality_forward, quality_reverse,
                       file_forward, file_reverse,
                       edit_func, iter_len, edit_len=1):
        """
        Driver method for creating edited FASTQ files.
        Args:
            edit_func: Function that performs the sequence edit (deletion, swap, or insertion)
            output_prefix: Prefix for output filenames
            edit_len: Length of edit (default=1)
        """
        # Determine output files based on unidirectional setting
        output_files = []
        if self.unidirectional:
            output_files.append(open(file_forward, 'w'))
        else:
            output_files.append(open(file_forward, 'w'))
            output_files.append(open(file_reverse, 'w'))

        for index in range(iter_len):
            # Generate headers
            header_forward = self.header_str.format(index=index, dir=1)
            header_reverse = self.header_str.format(index=index, dir=2)

            # Apply the edit
            edited_sequence = edit_func(sequence, index, edit_len, self._random)

            # Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
            sequence_forward = edited_sequence[:len(quality_forward)]
            sequence_reverse = edited_sequence[-len(quality_reverse):][::-1]

            # Write entries
            write_fastq_entry(output_files[0], header_forward,
                                    sequence_forward, quality_forward)

            if not self.unidirectional:
                write_fastq_entry(output_files[1], header_reverse,
                                        dna_complement(sequence_reverse),
                                        quality_reverse)
        # Ensure files are closed
        for f in output_files:
            f.close()

# Replaces the base pair at index 'index' with a random base pair
def single_pair_swap(sequence, index, r: random.Random):
    complement = dna_complement(sequence[index])
    list_sequence = list(sequence)
    list_sequence[index] = complement
    return ''.join(list_sequence)

# Inserts the complement of the pair at index 'index' after the pair
def single_pair_insertion(sequence, index, r: random.Random):
    complement = dna_complement(sequence[index])
    list_sequence = list(sequence)
    list_sequence.insert(index + 1, complement)
    return ''.join(list_sequence)

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
