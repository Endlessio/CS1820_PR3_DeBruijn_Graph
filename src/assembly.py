import enum
from contamination import Contamination
from correction import Correction
from preprocessing import get_reads, get_score_matrix, get_vector
import sys

'''
global variable
'''
ROOT_DIR = "./test_cases/assembly/"


class Assembly:
    def __init__(self):
        pass


def main(argv):
    # preprocess data
    vector = get_vector(ROOT_DIR+argv[1])
    reads = get_reads(ROOT_DIR+argv[2])
    contamination_k = int(argv[3])
    correction_k = int(argv[4])
    correction_d = int(argv[5])
    correction_t = int(argv[6])
    correction_mode = argv[7]
    matrix = get_score_matrix(ROOT_DIR+"unitary.m")

    # contamination to delete the reads
    contam = Contamination(reads, vector, contamination_k)
    kmer2idx_dict = contam.vector2kmer()
    _, res_idx = contam.end_match(kmer2idx_dict)
    
    res_idx = set(res_idx)
    contam_del_reads = []
    for idx, ele in enumerate(reads):
        if idx not in res_idx:
            contam_del_reads.append(ele)
    print(len(contam_del_reads), len(res_idx))
    assert len(contam_del_reads) == len(reads)-len(res_idx)


    # deleted reads for correction

    # corrected reads for debruijn

    pass

if __name__ == "__main__":
    main(sys.argv)