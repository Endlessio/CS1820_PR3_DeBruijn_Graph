from matplotlib import collections
from local_alignment import LocalAlignment
from debruijn import DeBruijn
from contamination import Contamination
from correction import Correction
from preprocessing import get_reads, get_score_matrix, get_vector, form_kmers, get_true_read
import sys
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import gc


'''
global variable
'''
ROOT_DIR = "./test_cases/assembly/"

def ablation_test(vector, reads):
    def contamination_k(k_range):
        res_reads = collections.defaultdict(list)
        res_length = []
        for k in k_range:
            contam = Contamination(reads, vector, k)
            kmer2idx_dict = contam.vector2kmer()
            _, contam_res_idx = contam.end_match(kmer2idx_dict)
            
            contam_res_idx = set(contam_res_idx)
            contam_del_reads = []
            for idx, ele in enumerate(reads):
                if idx not in contam_res_idx:
                    contam_del_reads.append(ele)
            res_length.append(len(contam_del_reads))
            res_reads[k] = contam_del_reads
        plt.plot(k_range, res_length)
        plt.xlabel("contamination k")
        plt.ylabel("reads amount after contamination")
        for x,y in zip(k_range, res_length):
            plt.text(x, y, y, ha="center", va="bottom", fontsize=10)
        plt.savefig('../res/ablation/contamination.png')
        
        return res_length, res_reads
    
    def correction(k_range, d_range, t_range, contamination_read):
        res_length = np.zeros(len(k_range), len(d_range), len(t_range))
        res = collections.defaultdict(list)
        for idx_k, k in enumerate(k_range):
            for idx_d, d in enumerate(d_range):
                for idx_t, t in enumerate(t_range):
                    correction = Correction(contamination_read, k, t, d)
                    kmer_list, tot_kmer_dict = correction.form_kmer(k)
                    infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, t)
                    closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, d, k)

                    # choose your algorithm
                    corrected_res, _ = correction.stack_replace(closest_pair_dict)
                    res_length[idx_k][idx_d][idx_t] = len(corrected_res)
                    res[(idx_k, idx_d, idx_t)] = corrected_res
                    # print("corrected read number:", len(corrected_res))

        # formulate k, d, t, v
        k_list = []
        d_list = []
        t_list = []
        v_list = []
        for key, val in res.items():
            k_list.append(key[0])
            d_list.append(key[1])
            t_list.append(key[2])
            v_list.append(len(val))
        '''
        fig = plt.figure()   
        ax = fig.add_subplot(projection='3d')
        ax.scatter(k_list, d_list, v_list)
        ax.set_xlabel('k')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        '''




    contam_res_length, contam_res_reads = contamination_k([2,4,6,8,10,12,14,16])
    print("after contamination, the reads amount:", contam_res_length)





def main(argv):
    # preprocess data
    vector = get_vector(ROOT_DIR+"vector.txt")
    print(f'vector length, {len(vector)}')
    reads = get_reads(ROOT_DIR+"sampleReads.txt")
    print(f'reads number, {len(reads)}')
    
    '''
    # Uncomment this part if you want to re-run the ablation test
    # or the results can be found in res/ablation/contamination.png
    # --------------------------------------------------------------------------
    ablation_test(vector, reads)
    '''

    contamination_k = int(argv[3]) #  超过k的连续subsequence被认为是contamination
    correction_k = int(argv[4]) # k-mer
    correction_d = int(argv[5]) # 替换infrequent时，最多差了d个position
    correction_t = int(argv[6]) # 出现次数小于t的k-mer被认为是infrequent
    graph_k = int(argv[7]) # k-mer length when construct the graph
    correction_mode = argv[8]
    matrix = get_score_matrix(ROOT_DIR+"unitary.m")

    # contamination to delete the reads
    contam = Contamination(reads, vector, contamination_k)
    kmer2idx_dict = contam.vector2kmer()
    _, contam_res_idx = contam.end_match(kmer2idx_dict)
    
    contam_res_idx = set(contam_res_idx)
    contam_del_reads = []
    for idx, ele in enumerate(reads):
        if idx not in contam_res_idx:
            contam_del_reads.append(ele)
    # print("after contamination, left read amount,", len(contam_del_reads))
    assert len(contam_del_reads) == len(reads)-len(contam_res_idx)


    # deleted reads for correction
    correction = Correction(contam_del_reads, correction_k, correction_t, correction_d)

    kmer_list, tot_kmer_dict = correction.form_kmer(correction_k)
    infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, correction_t)
    closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, correction_d, correction_k)

    # choose your algorithm
    if correction_mode == "stack":
        corrected_res, _ = correction.stack_replace(closest_pair_dict)
    elif correction_mode == "simple":
        corrected_res, _ = correction.simple_replace(closest_pair_dict)
    elif correction_mode == "naive":
        corrected_res, _ = correction.naive_replace(closest_pair_dict)
    elif correction_mode == "merge":
        corrected_res, _ = correction.opt_merge_replace(closest_pair_dict, infrequent_kmer_dict)
    else:
        print("wrong correction mode")
        quit(0)
    print("corrected read number:", len(corrected_res))

    # get divide reads to k-mer
    graph_reads = list(set(form_kmers(corrected_res, graph_k)))
    print("check k-mer number before and after", len(corrected_res), len(graph_reads))

    del correction
    del corrected_res
    gc.collect()

    # corrected reads for debruijn
    my_graph = DeBruijn(graph_reads)
    my_graph.construct_graph()
    # print("done construction")
    my_graph.merge_singleton()
    # print("done merge")
    start_nodes = my_graph.get_start_node()
    # print("done get start")
    res, _ = my_graph.find_path(start_nodes)
    # print("done get path")
    print("number of isolates:", nx.number_of_isolates(my_graph.graph), "number of nodes:", my_graph.graph.number_of_nodes())
    my_graph.print_out(res)



if __name__ == "__main__":
    main(sys.argv)