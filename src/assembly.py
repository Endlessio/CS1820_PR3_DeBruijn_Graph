import contamination
import correction
from preprocessing import get_reads, get_score_matrix, get_vector
import sys

'''
global variable
'''
ROOT_DIR = "../test_cases/assembly/"


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

    # contamination to delete the reads

    # deleted reads for correction

    # corrected reads for debruijn

    pass

if __name__ == "__main__":
    main(sys.argv)