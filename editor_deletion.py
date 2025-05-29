import quality_score_approximation
import sys

DELETION_COUNT = 1 if len(sys.argv) < 1 else sys.argv[1]

def main():
    # Read in the file containing reference sequence
    reference_file = "sperm.referenceseq.fa"
    refence_sequence = open_file_read(reference_file)
    #Get rid of the first line so that only raw sequence data remains
    refence_sequence = refence_sequence.splitlines()[1]
    print(f'reference sequence:\n{refence_sequence}\n')
    
    #Store the file path of this original reads
    reads_file_path_1 = 'reference_reads_compressed/read1.fastq.gz'
    reads_file_path_2 = 'reference_reads_compressed/read2.fastq.gz'
    
    #Get the average quality scores for each file (sampled)
    qual_string_1 = quality_score_approximation.get_average_quality(reads_file_path_1)
    qual_string_2 = quality_score_approximation.get_average_quality(reads_file_path_2)
    print(f'Average quality for forward read:\n{qual_string_1}')
    print(f'Average quality for backward read:\n{qual_string_2}')
        
    #Create the output files
    create_simulated_fastqs(refence_sequence, qual_string_1, qual_string_2)
    
# Creates fastq files with deletions for each base pair
def create_simulated_fastqs_deletion(sequence, quality_forward, quality_backward):
    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
    with open('single_deletion_forward.fastq', "w") as f1:
        with open('single_deletion_reverse.fastq', "w") as f2:
            for index in range(len(sequence) - DELETION_COUNT):
                header_forward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 1:N:0:GATCAGAT+TCCGCGAA'
                header_backward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 2:N:0:GATCAGAT+TCCGCGAA'
                edited_sequence = sequence[:index] + sequence[index + DELETION_COUNT:]
                #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                forward_sequence = edited_sequence[:len(quality_forward)]
                backward_sequence = edited_sequence[-len(quality_backward):][::-1]
                
                #Write the contents to the appropriate files
                f1.write(f'{header_forward}\n{forward_sequence}\n+\n{quality_forward}\n')
                f2.write(f'{header_backward}\n{dna_complement(backward_sequence)}\n+\n{quality_backward}\n')
                
# Creates fastq files with replacements for each base pair
def create_simulated_fastqs_replacement(sequence, quality_forward, quality_backward):
    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
    with open('single_swap_forward.fastq', "w") as f1:
        with open('single_swap_reverse.fastq', "w") as f2:
            for index in range(len(sequence) - 1):
                header_forward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 1:N:0:GATCAGAT+TCCGCGAA'
                header_backward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 2:N:0:GATCAGAT+TCCGCGAA'
                edited_sequence = single_pair_swap(sequence, index)
                #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                forward_sequence = edited_sequence[:len(quality_forward)]
                backward_sequence = edited_sequence[-len(quality_backward):][::-1]
                
                #Write the contents to the appropriate files
                f1.write(f'{header_forward}\n{forward_sequence}\n+\n{quality_forward}\n')
                f2.write(f'{header_backward}\n{dna_complement(backward_sequence)}\n+\n{quality_backward}\n')

# Creates fastq files with insertions between each pair of base pairs
def create_simulated_fastqs_insertion(sequence, quality_forward, quality_backward):
    #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
    with open('single_insertion_forward.fastq', "w") as f1:
        with open('single_insertion_reverse.fastq', "w") as f2:
            for index in range(len(sequence) - 1):
                header_forward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 1:N:0:GATCAGAT+TCCGCGAA'
                header_backward = f'@M05774:236:000000000-LKC5Y:1:1107:27241:{str(index)} 2:N:0:GATCAGAT+TCCGCGAA'
                edited_sequence = single_pair_insertion(sequence, index)
                #Create forward and backward sequences with slicing. Backward sequence is reversed with slicing.
                forward_sequence = edited_sequence[:len(quality_forward)]
                backward_sequence = edited_sequence[-len(quality_backward):][::-1]
                
                #Write the contents to the appropriate files
                f1.write(f'{header_forward}\n{forward_sequence}\n+\n{quality_forward}\n')
                f2.write(f'{header_backward}\n{dna_complement(backward_sequence)}\n+\n{quality_backward}\n')

# Replaces
def single_pair_swap(sequence, index):
    complement = dna_complement(sequence[index])
    list_sequence = list(sequence)
    list_sequence[index] = complement
    return ''.join(list_sequence)

# Inserts the complement of the pair at index 'index' after the pair
def single_pair_insertion(sequence, index):
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

def open_file_read(path):
    try:
        with open(path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print(f"Error: The file '{path}' was not found.")
        exit(1)

if __name__ == "__main__":
    main()
