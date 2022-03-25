from collections import defaultdict


def get_reads(file):
    reads = []
    with open(file) as f:
        for line in f.readlines():
            cur = line.strip()
            reads.append(cur)
    return reads


def get_vector(file):
    with open(file) as f:
        vector = f.readline().strip()
    return vector


def get_score_matrix(file):
    matrix = defaultdict(dict)
    with open(file) as f:
        chars = f.readline().strip().split()[1:]
        temp = []
        for _ in range(len(chars)):
            temp.append(f.readline()[1:].strip().split())
        row = len(temp)
        col = len(temp[0])
        for i in range(row):
            for j in range(col):
                matrix[chars[i]][chars[j]] = int(temp[i][j])

    return matrix


