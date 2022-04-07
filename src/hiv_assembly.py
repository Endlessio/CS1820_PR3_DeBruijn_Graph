import os
from threading import local
from local_alignment import LocalAlignment, global_alignment
from debruijn import DeBruijn
from contamination import Contamination
from correction import Correction
from preprocessing import get_reads, get_score_matrix, get_vector, form_kmers, get_true_read
import matplotlib.pyplot as plt
import random
from collections import defaultdict
import heapq
from matplotlib import cm
import numpy as np
import pickle


'''
global variables
'''
ROOT_DIR = "../test_cases/assembly/"
HIV_GENOME_PATH = "../test_cases/assembly/hiv-1_genome.txt"
HIV_READS_PATH = "../test_cases/assembly/hiv-1_reads.txt"
VECTOR_PATH = "../test_cases/assembly/vector.txt"
MATRIX_PATH = "../test_cases/assembly/unitary.m"

CONTAMINATION_k = 6
CORRECTION_k = 17
CORRECTION_t = 4
CORRECTION_d = 2
GRAPH_k = [10,20,30,40,50]
COVERAGE = [4,6,8,10]
MATRIX = get_score_matrix(MATRIX_PATH)


def generate_figures():
    a_range = []
    k_range = []

    avg_length = []
    global_align_score = []
    local_align_score = []

    zip_form = []

    file_list = os.listdir("../res/contigs")

    for file_name in file_list:
        a = int(file_name.split("_")[2])
        k = int(file_name.split("_")[-1][:-4])
        with open("../res/alignment_analysis/"+file_name[:-4]+".pkl", 'rb') as f:
            data = pickle.load(f)
            length = data[1]
            g_score = data[2]
            l_score = data[3]

            a_range.append(a)
            k_range.append(k)
            avg_length.append(int(length))
            global_align_score.append(g_score)
            local_align_score.append(l_score)
            zip_form.append((a,k,length,g_score, l_score))
        f.close()
    print(a_range)
    print(k_range)
    print(avg_length)
    print(global_align_score)
    print(local_align_score)
    print(sorted(zip_form))
    
    # # Creating figure
    fig1 = plt.figure()
    ax1 = plt.axes(projection ="3d")
    

    ax1.scatter3D(a_range, k_range, avg_length, cmap=cm.jet)
    plt.title("coverage, graph k-mer, and average contig length")
    ax1.set_xlabel('Coverage Range')
    ax1.set_ylabel('Graph K range')
    ax1.set_zlabel('Average Contig Length')
    plt.savefig("../res/a_k_avg-contig-len.png")


    # ax1.plot_trisurf(a_range, k_range, global_align_score)
    # plt.title("coverage, graph k-mer, and avg global alignment score")
    # ax1.set_xlabel('Coverage Range')
    # ax1.set_ylabel('Graph K range')
    # ax1.set_zlabel('Global Alignment Score')
    # plt.savefig("../res/a_k_avg-global.png")


    # ax1.plot_trisurf(a_range, k_range, local_align_score, cmap=cm.jet)
    # plt.title("coverage, graph k-mer, and avg local alignment score")
    # ax1.set_xlabel('Coverage Range')
    # ax1.set_ylabel('Graph K range')
    # ax1.set_zlabel('Local Alignment Score')
    # plt.savefig("../res/a_k_avg-local.png")



def sample(reads, number):
    return random.choices(reads, k=number)

def contig_number(N, L):
    return N*np.exp(-(N*L/9181))

def calculate_sample(L, G, coverage):
    return coverage*G/L

def alignment_analysis(contig_file_path):
    def get_avg_length(contigs_path):
        length_sum = 0
        cnt = 0
        contigs = []
        concatenate = ''
        with open(contigs_path, "r") as f:
            for line in f.readlines():
                line = line.strip()
                if line[0] != ">":
                    if length_sum < 10000:
                        contigs.append(line)
                        concatenate += line
                        length_sum += len(line)
                        cnt += 1
                    else:
                        break
        f.close()
        return length_sum/cnt, contigs, concatenate

    true_genome = get_true_read(HIV_GENOME_PATH).strip()
    contigs_list = os.listdir(contig_file_path)

    for contigs_name in contigs_list:
        single_align_score_list = []
        average_length, contigs, concatenate = get_avg_length("../res/contigs/"+contigs_name)
        global_align_score = global_alignment(concatenate, true_genome)
        norm_global_align_score = global_align_score/9181
        for contig in contigs:
            align = LocalAlignment([contig, true_genome], MATRIX, -1)
            align.local_alignment()
            single_align_score_list.append(align.imax)
        to_write = [single_align_score_list, average_length, global_align_score, norm_global_align_score, contigs]
        print(f"{contigs_name[:-4]}: {to_write}")
        p = open("../res/alignment_analysis/"+contigs_name[:-4]+".pkl","wb")
        # write the python object (dict) to pickle file
        pickle.dump(to_write, p)
        p.close()


    '''
    alignment_score = []
    alignment_length = []
    for_global_align = ""

    # alignment
    for_global_align += contig[1]
    align = LocalAlignment([contig[1], true_genome], MATRIX, -1)
    align.local_alignment()
    alignment_length.append(len(contig[1]))
    alignment_score.append(align.imax/len(contig[1]))
    # print(alignment_score, sum(alignment_score)/len(alignment_score))

    global_align[(a, cur_k)] = [global_alignment(for_global_align, true_genome), sum(alignment_score)/len(alignment_score), alignment_score, alignment_length]
    print(global_align)

    # create a binary pickle file 
    p = open("../test_cases/assembly/global_alignment_final.pkl","wb")
    # write the python object (dict) to pickle file
    pickle.dump(global_align,p)
    # close file
    p.close()
    '''
    

def ablation_test(vector, reads):
    def contamination_k(k_range):
        res_reads = defaultdict(list)
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
        plt.savefig('../res/ablation/hiv_contamination.png')
        plt.close()
        
        return res_length, res_reads

    def coverage(coverage_list):
        L = 50
        G = 9181
        N = []
        for c in coverage_list:
            N.append(c*G//L)

        plt.plot(coverage_list, N)
        plt.xlabel("coverage")
        plt.ylabel("reads amount need to sample")
        for x,y in zip(coverage_list, N):
            plt.text(x, y, y, ha="center", va="bottom", fontsize=10)
        plt.savefig('../res/ablation/hiv_coverage.png')
        plt.close()
        return N, coverage_list

    contam_res_length, contam_res_reads = contamination_k([2,4,6,8,10,12,14,16])
    print("after contamination, the reads amount:", contam_res_length)

    N_list, coverage_list = coverage([2,4,6,8,10,12])
    print("for coverage:", coverage_list, ", the N is:", N_list)

def main():
    # get true genome
    true_genome = get_true_read(HIV_GENOME_PATH)
    print("the true reads length is:", len(true_genome))
    if not os.path.exists("../res/true_reads/hiv_true_genome.txt"):
        with open("../res/true_reads/hiv_true_genome.txt", "w") as f:
            f.write("> true reads")
            f.write("\n")
            f.write(true_genome)

    # get reads
    reads = get_reads(HIV_READS_PATH)
    print("each reads length is:", len(reads[0]), "; reads amount is:", len(reads))

    # get vector
    vector = get_vector(VECTOR_PATH)

    '''
    # Uncomment this part if you want to redo the ablation test, contamination and correction
    # ---------------------------------------------------------------------------------------
    # ablation test
    # result:
        # contamination k: 10% contamination, 9681-918 = 8236; choose 6, get 8763
        # from project 2, we have the local optimal hyper-parameter
            # CORRECTION_k = 17
            # CORRECTION_t = 4
            # CORRECTION_d = 2
        # coverage: calculate coverage of corresponding number
            # [2, 4, 6, 8, 10, 12] -> [367, 734, 1101, 1468, 1836, 2203]
        # get multiple res with certain threshold
    ablation_test(vector, reads)


    # apply contamination
    contam = Contamination(reads, vector, CONTAMINATION_k)
    kmer2idx_dict = contam.vector2kmer()
    _, contam_res_idx = contam.end_match(kmer2idx_dict)
    
    contam_res_idx = set(contam_res_idx)
    contam_del_reads = []
    for idx, ele in enumerate(reads):
        if idx not in contam_res_idx:
            contam_del_reads.append(ele)
    print("after contamination read amount", len(contam_del_reads))
    assert len(contam_del_reads) == len(reads)-len(contam_res_idx)


    # apply correction
    correction = Correction(contam_del_reads, CORRECTION_k, CORRECTION_t, CORRECTION_d)

    kmer_list, tot_kmer_dict = correction.form_kmer(CORRECTION_k)
    infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, CORRECTION_t)
    closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, CORRECTION_d, CORRECTION_k)
    corrected_res, _ = correction.stack_replace(closest_pair_dict)
    print("corrected read number:", len(corrected_res))

    # save corrected reads to txt
    with open('../test_cases/assembly/hiv_corrected_reads.txt', 'w') as f:
        for r in corrected_res:
            f.write(r)
            f.write('\n')
    '''



    '''
    # Uncomment this part if you want to redo the De Bruijn genome assembly process
    # or the results are stored in res/contigs/
    # ---------------------------------------------------------------------------------------
    # sample
    reads = get_reads("../test_cases/assembly/hiv_corrected_reads.txt")

    for a in COVERAGE:
        for cur_k in GRAPH_k:
            sampled_reads = []
            sampled_amount = int(calculate_sample(cur_k, 9181, a))
            sampled_reads = sample(reads, sampled_amount)
            # get divide reads to k-mer
            graph_reads = list(set(form_kmers(sampled_reads, cur_k)))
            print("check k-mer number before and after", cur_k, len(sampled_reads), len(graph_reads))

            # deBruijn graph assembly
            graph = DeBruijn(graph_reads)
            graph.construct_graph()
            graph.merge_singleton()
            start_nodes = graph.get_start_node()
            # graph.detect_cycle()
            max_res, res = graph.find_path(start_nodes)
            threshold = int(contig_number(len(graph_reads), cur_k))
            contigs = heapq.nlargest(int(len(res)*0.07), res)
            print(len(contigs))

            print("after genome assembly by deBruijn:", [len(i[1]) for i in contigs])


            with open(f"../res/contigs/hiv_contig_{a}_{cur_k}.txt", "w") as f:
                for idx,contig in enumerate(contigs):
                    # write
                    f.write(">"+str(idx)+"\n")
                    f.write(contig[1])
                    f.write("\n")
    '''          

    '''
    # Uncomment this if you want to re-run the alignment analysis result
    # or the results are stored in res/alignment_analysis/
    # ---------------------------------------------------------------------------------------
    alignment_analysis("../res/contigs/")
    '''

    '''
    # Uncomment this if you want to re-run the analysis figure generation again
    # or the results are stored in res/
    # ---------------------------------------------------------------------------------------
    '''
    generate_figures()




main()