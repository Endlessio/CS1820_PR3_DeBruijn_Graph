def get_reads(reads_file):
    reads = []
    with open(reads_file) as f:
        for line in f.readlines():
            cur = line.strip()
            reads.append(cur)
    return reads


