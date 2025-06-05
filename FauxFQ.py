import pathlib
import quality_score_approximation
import argparse
from fq_editor import Fqeditor

parser = argparse.ArgumentParser(
                    prog='FauxFastQ',
                    description='Creates a fastQ file with simulated edits at each point in the reference sequence.',
                    epilog='https://github.com/DrewShapiro5/FauxFastQ')

parser.add_argument('reference_fa', type=pathlib.Path,)
parser.add_argument('output_forward', type=pathlib.Path, help="Path to simulated output file (forward sequence)")
parser.add_argument('output_reverse', type=pathlib.Path, help="Path to simulated output file (reverse sequence)")
parser.add_argument('-n', '--edit_length', type=int, default=1, help='Amount of base pairs in the edits.')
parser.add_argument('--unidirectional', action='store_true', help="Use unidirectional sequencing")
parser.add_argument('--reference_forward', type=pathlib.Path, help="Path to forward reference fastq file (for quality score sampling)")
parser.add_argument('--reference_reverse', type=pathlib.Path, help="Path to reverse reference fastq file (for quality score sampling)")
parser.add_argument('-s', '--seed', type=int, help="Seed for random base generator")
parser.add_argument('-d', '--deletion', action='store_true', help="Deletion mode")
parser.add_argument('-i', '--insertion', action='store_true', help="Insertion mode")
parser.add_argument('-r', '--replacement', action='store_true', help="Replacement mode")
parser.add_argument('--append', action='store_true', help="For use with multiple editing modes. Creates a fasta/q file with edits of each specified mode.")
parser.add_argument('--sample_size', type=int, default=500, help="Sample size for quality sampling from reference files")
parser.add_argument('--use_quality', type=bool, default=True, help="Output files will contain quality scores by default (fastq format). If set to false, output files will be in fasta format.")
parser.add_argument('-l', '--read_length', type=int, required=False, help="Use quality score approximation")

args = parser.parse_args()
print(args)

def main():
    reference_file = args.reference_fa
    # Read in the file containing reference sequence
    # Get rid of the first line so that only raw sequence data remains
    reference_sequence = open_file_read(reference_file).splitlines()[1]
    print(f'reference sequence:\n{reference_sequence}\n')

    # Initialize variables
    seq_hash = hash(reference_sequence)
    seed = args.seed if args.seed else seq_hash
    unidirectional = args.unidirectional if args.unidirectional else False
    file_mode = 'a'

    # Check if output files are the same
    if args.output_forward == args.output_reverse and not unidirectional:
        print(f'Output file paths are the same: {args.output_forward}\nNo outputs will be created.')
        exit(1)

    # Get the average quality scores for each reference file (sampled)
    qual_string_forward = None
    qual_string_reverse = None
    if args.reference_forward:
        reference_forward = args.reference_forward
        qual_string_forward = quality_score_approximation.get_average_quality(str(reference_forward), args.sample_size)
        print(f'Average quality for forward read:\n{qual_string_forward}\n')
    elif args.read_length:
        qual_string_forward = max_quality_string(args.read_length)
    else:
        print('No forward reference file provided for quality sampling.'
              'Provide a reference fastq file with --reference_forward <file>, or specify read length with -l <length>')
        exit(1)

    if args.reference_reverse and not unidirectional:
        reference_reverse = args.reference_reverse
        qual_string_reverse = quality_score_approximation.get_average_quality(str(reference_reverse), args.sample_size)
        print(f'Average quality for reverse read:\n{qual_string_reverse}\n')
    elif args.read_length and not unidirectional:
        qual_string_reverse = max_quality_string(args.read_length)
    elif not unidirectional:
        print('No reverse reference file provided for quality sampling.'
              'Provide a reference fastq file with --reference_reverse <file>, or specify read length with -l <length>')
        exit(1)

    # If not in append mode, delete the previous contents of the output files
    if not args.append:
        with open(args.output_forward, 'w') as _:
            pass
        if not args.unidirectional:
            with open(args.output_reverse, 'w') as _:
                pass

    total_edit_modes = 0
    for mode in [args.deletion, args.insertion, args.replacement]:
        total_edit_modes += 1 if mode else 0
    if total_edit_modes < 1 and not args.append:
        print('No edit mode specified. Exiting.')
        exit(0)

    # Create the fastq editor
    fqe = Fqeditor(seed=seed, unidirectional=unidirectional, file_mode=file_mode, file_forward=str(args.output_forward), file_reverse=str(args.output_reverse))
    # TODO: runs for fasta, figure out how to index headers correctly with append mode (use sequence hash + edit parameters as the machine id?)
    if args.deletion:
        print('Performing deletions...\n')
        fqe.create_simulated_fastqs_deletion(sequence=reference_sequence,
                                             quality_forward=qual_string_forward,
                                             quality_reverse=qual_string_reverse,
                                             edit_len=args.edit_length,)
    if args.insertion:
        print('Performing insertions...\n')
        fqe.create_simulated_fastqs_insertion(sequence=reference_sequence,
                                              quality_forward=qual_string_forward,
                                              quality_reverse=qual_string_reverse,
                                              edit_len=args.edit_length,)
    if args.replacement:
        print('Performing replacements...\n')
        fqe.create_simulated_fastqs_replacement(sequence=reference_sequence,
                                                quality_forward=qual_string_forward,
                                                quality_reverse=qual_string_reverse,
                                                edit_len=args.edit_length,)

def open_file_read(path):
    try:
        with open(path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print(f"Error: The file '{path}' was not found.")
        exit(1)

def max_quality_string(length):
    return '!' * length

if __name__ == "__main__":
    main()