import pathlib
import fq_editor
import quality_score_approximation
import argparse
import pyfastx

from fq_editor import Fqeditor

parser = argparse.ArgumentParser(
                    prog='FauxFastQ',
                    description='Creates a fastQ file with simulated edits at each point in the reference sequence.',
                    epilog='https://github.com/DrewShapiro5/FauxFastQ')

parser.add_argument('reference_fa', type=pathlib.Path,)
parser.add_argument('-n', '--edit_length', type=int, default=1, help='Amount of base pairs in the edits.')
parser.add_argument('--unidirectional', action='store_true', help="Use unidirectional sequencing")
parser.add_argument('--reference_forward', type=pathlib.Path, help="Path to forward reference fastq file (for quality score sampling)")
parser.add_argument('--reference_reverse', type=pathlib.Path, help="Path to reverse reference fastq file (for quality score sampling)")
parser.add_argument('--output_forward', type=pathlib.Path, help="Path to simulated output file (forward sequence)")
parser.add_argument('--output_reverse', type=pathlib.Path, help="Path to simulated output file (reverse sequence)")
parser.add_argument('-s', '--seed', type=int, help="Seed for random number generator")
parser.add_argument('-d', '--deletion', action='store_true', help="Deletion mode")
parser.add_argument('-i', '--insertion', action='store_true', help="Insertion mode")
parser.add_argument('-r', '--replacement', action='store_true', help="Replacement mode")
parser.add_argument('--append', action='store_true', help="For use with multiple editing modes. Creates a fasta/q file with edits of each specified mode.")

args = parser.parse_args()
print(args)

def main():
    reference_file = args.reference_fa
    # Read in the file containing reference sequence
    # Get rid of the first line so that only raw sequence data remains
    reference_sequence = open_file_read(reference_file).splitlines()[1]
    print(f'reference sequence:\n{reference_sequence}\n')

    seq_hash = hash(reference_sequence)
    seed = args.seed if args.seed else seq_hash
    unidirectional = args.unidirectional if args.unidirectional else False

    qual_string_forward = ''
    if args.reference_forward:
        # Store the file path of this original reads
        reference_forward = args.reference_forward
        qual_string_forward = quality_score_approximation.get_average_quality(reference_forward)
        print(f'Average quality for forward read:\n{qual_string_forward}')

    qual_string_reverse = ''
    if args.reference_reverse:
        reference_reverse = args.reference_reverse
        # Get the average quality scores for each file (sampled)
        qual_string_reverse = quality_score_approximation.get_average_quality(reference_reverse)
        print(f'Average quality for reverse read:\n{qual_string_reverse}')

    file_mode = 'a' if args.append else 'w'
    fqe = Fqeditor(seed=seed, unidirectional=unidirectional, file_mode=file_mode)
    total_edit_modes = 0
    if args.deletion:
        # TODO: for runs without quality score, specify read length
        fqe.create_simulated_fastqs_deletion(sequence=reference_sequence,
                                             quality_forward=qual_string_forward,
                                             quality_reverse=qual_string_reverse,
                                             edit_len=args.edit_length,
                                             file_forward=args.output_forward,
                                             file_reverse=args.output_reverse)
        total_edit_modes += 1
    if args.insertion:
        # TODO: for runs without quality score, specify read length
        fqe.create_simulated_fastqs_insertion(sequence=reference_sequence,
                                              quality_forward=qual_string_forward,
                                              quality_reverse=qual_string_reverse,
                                              edit_len=args.edit_length,
                                              file_forward=args.output_forward,
                                              file_reverse=args.output_reverse)
        total_edit_modes += 1
    if args.replacement:
        # TODO: for runs without quality score, specify read length
        fqe.create_simulated_fastqs_replacement(sequence=reference_sequence,
                                                quality_forward=qual_string_forward,
                                                quality_reverse=qual_string_reverse,
                                                edit_len=args.edit_length,
                                                file_forward=args.output_forward,
                                                file_reverse=args.output_reverse)
        total_edit_modes += 1
    if total_edit_modes > 1 and not args.append:
        print('Warning: multiple edit modes specified without appending. Use --append or only one edit mode.')



def open_file_read(path):
    try:
        with open(path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print(f"Error: The file '{path}' was not found.")
        exit(1)

if __name__ == "__main__":
    main()