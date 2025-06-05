import pyfastx
import random
import sys

def get_average_quality(file_path, sample_size):
    fa = pyfastx.Fastq(file_path)
    reads_count = len(fa)
    print('Average quality score outputs:')
    
    print(f'reads: {reads_count}')

    #These should all be the same
    print(f'mean length: {fa.avglen}')
    print(f'min length: {fa.minlen}')
    print(f'max length: {fa.maxlen}')
    read_length = fa.minlen
    sample_size = min(sample_size, reads_count)
    sample_indices = []

    #Create a list of random indices
    for i in range(sample_size):
        random_index = random.randint(0, reads_count)
        while random_index in sample_indices:
            random_index = random.randint(0, reads_count)
        sample_indices.append(random_index)

    #This array will store the sum of all quality scores
    quality_scores_sum = [0] * read_length

    #Add the quality of all reads in the sample
    for index in sample_indices:
        r = fa[index]
        quality_score = r.quali
        for i in range(read_length):
            quality_scores_sum[i] += quality_score[i]

    #Normalize the quality scores
    quality_scores_mean = []
    for i in range(read_length):
        quality_scores_mean.append(quality_scores_sum[i] / sample_size)
        
    #Create a string representing the quality score true to fastq format
    #The following string is found on the fastq wikipedia page
    quality_chars ='!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

    quality_score_string = ''
    for val in quality_scores_mean:
        quality_score_string += quality_chars[round(val)]

    return quality_score_string