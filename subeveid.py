#!/home/max/bin/anaconda3/bin/python

import sys, os
from collections import OrderedDict
import pickle
import numpy as np


############## CLI ARGS #######################
min_s = int(sys.argv[1])            # min size to consider piRNA
max_s = int(sys.argv[2])            # max size to consider piRNA
outpath = "/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/subeveid_out"     # path to write to
w_size = 500            # chunk size
# load sRNA BED file 

# sliding windows - length, stepsize
##################### GLOBALS #################################
EVEs_dict = {} # {name:[chr, startpos, stoppos, length, strand]}

### results data structure #####
# {EVE_name : [chunk1 : [startpos, stoppos, piRNA_count_p, max_depth_p, mean_depth_p, sum_depth_p, {pos1+:count, pos2+:count....pos_n+:count}, piRNA_count_m, max_depth_m, mean_depth_m, sum_depth_m, {pos1-:count, pos2-:count....pos_n-:count}], chunk2 : [...]]}
piRNA_EVE_p = {} # collect piRNAs per EVE mapping in pos orientation
piRNA_EVE_m = {} # collect piRNAs per EVE mapping in neg orientation

### chromosome dict structure ###
# {startpos: [[read1],[read2],...,[readn]]}
# [read1] = name, stoppos, orientation
chr1_reads_p = {}
chr2_reads_p = {}
chr3_reads_p = {}
chr1_reads_m = {}
chr2_reads_m = {}
chr3_reads_m = {}



###################### FUNCTION ################################

def load_read(read_entry, chrom, start, stop):   
    data = chrom.get(startpos)
    if data == None:
        chrom.update({startpos:[[read_entry[3], stoppos, read_entry[5]]]})
    else:
        chrom[startpos].append([read_entry[3], stoppos, read_entry[5]])                   
        # print(f'{startpos}:\t{chr1_reads_p.get(startpos)}')
def save_pickles():
    global chr1_reads_p, chr2_reads_p, chr3_reads_p, chr1_reads_m, chr2_reads_m, chr3_reads_m
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr1_reads_p.pkl', 'wb') as chr1p:
        pickle.dump(chr1_reads_p, chr1p, pickle.HIGHEST_PROTOCOL)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr2_reads_p.pkl', 'wb') as chr2p:
        pickle.dump(chr2_reads_p, chr2p, pickle.HIGHEST_PROTOCOL)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr3_reads_p.pkl', 'wb') as chr3p:
        pickle.dump(chr3_reads_p, chr3p, pickle.HIGHEST_PROTOCOL)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr1_reads_m.pkl', 'wb') as chr1m:
        pickle.dump(chr1_reads_m, chr1m, pickle.HIGHEST_PROTOCOL)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr2_reads_m.pkl', 'wb') as chr2m:
        pickle.dump(chr2_reads_m, chr2m, pickle.HIGHEST_PROTOCOL)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr3_reads_m.pkl', 'wb') as chr3m:
        pickle.dump(chr3_reads_m, chr3m, pickle.HIGHEST_PROTOCOL)

def count_pi_perpos(chrom, posit, chunk_lim):
    ''' Count piRNA depth per position and save info for current chunk'''
    global chunk_depth_p, pi_nums_p
    data = chrom.get(int(posit))
    readcount = 0
    if data != None:                          
        for entry in data:   
            # print(f'{read[0]}\t{read[1]}')
            if int(entry[1]) <= chunk_lim:    # only count reads that dont exceed the chunk border shall be accounted for
                readcount += 1
        if readcount > 0:                      # only if valid reads counted, increase individual piRNA count
            pi_nums_p += 1                          
    chunk_depth_p.update({int(posit):readcount})   # save read counts for pos in chunk

def count_pi_perneg(chrom, posit, chunk_lim):
    ''' Count piRNA depth per position and save info for current chunk'''
    global chunk_depth_m, pi_nums_m
    data = chrom.get(int(posit))
    readcount = 0
    if data != None:                          
        for entry in data:   
            # print(f'{entry[0]}\t{entry[1]}')
            if int(entry[1]) >= chunk_lim:    # only count reads that dont exceed the chunk border shall be accounted for
                readcount += 1
        if readcount > 0:                      # only if valid reads counted, increase individual piRNA count
            pi_nums_m += 1                              
    chunk_depth_m.update({int(posit):readcount})   # save read counts for pos in chunk

def calc_metrics(pi_rnas, chunk_depth):
    ''' Calculate metrics per chunk from chunk_depth file'''
    if pi_rnas != 0:                    # if piRNAs found on this chunk, calculate mean and max depth
        depths = list(chunk_depth.values())
        depths = np.array(depths)
        max_depth = np.amax(depths)
        sum_depth = np.sum(depths)
        mean_depth = sum_depth/pi_rnas  # depth per piRNA; how many reads cover one starting piRNA on average?
        # print(f'EVE\tchunk_name\tchunkstart\tchunkstop\tmean_depth\tmax_depth\tpi_nums')                    
    else:
        mean_depth = 0
        max_depth = 0   
        sum_depth = 0 
    # print(len(list(chunk_depth.values())))
    return sum_depth, mean_depth, max_depth 
###################### MAIN ####################################
# Load the EVE annotation
print("Loading EVEs...")
with open("/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/flavi_NIRVS_AeAe.csv", 'r') as EVEs:
    linecount = 0
    for line in EVEs:
        linecount += 1
        if linecount >= 3 and line.startswith("Flaviviridae"):
            entry = line.split("\t")
            EVEs_dict.update({entry[1]:[entry[2],entry[3],entry[4],entry[5],entry[6]]})
print("Loading complete.")

# Load BED file and establish data structure:
print("Loading piRNAs...")
if os.path.isfile(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr1_reads_p.pkl'): # ask for one pickle, imply all are there!
    # load the piRNAs from pickles, faster script running
    print("Load from pickles...")
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr1_reads_p.pkl', 'rb') as load_obj1:
        chr1_reads_p = pickle.load(load_obj1)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr2_reads_p.pkl', 'rb') as load_obj2:
        chr2_reads_p = pickle.load(load_obj2)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr3_reads_p.pkl', 'rb') as load_obj3:
        chr3_reads_p = pickle.load(load_obj3)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr1_reads_m.pkl', 'rb') as load_obj4:
        chr1_reads_m = pickle.load(load_obj4)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr2_reads_m.pkl', 'rb') as load_obj5:
        chr2_reads_m = pickle.load(load_obj5)
    with open(f'/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/chr3_reads_m.pkl', 'rb') as load_obj6:
        chr3_reads_m = pickle.load(load_obj6)
else:
    print("Load from bed file...")
    with open("/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/RDVJ106.bed", 'r') as reads: #"/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/flavi_EVEs/RDVJ106.bed"
        for line in reads:
            entry = line.rstrip("\n").split("\t")
            length = int(entry[2]) - int(entry[1])      # sRNA length
            if entry[5] == "+":
                startpos = int(entry[1])
                stoppos = int(entry[2])
                # print(f'{line}{length}')  
                if entry[0] == "1" and length >= min_s and length <= max_s:  
                    load_read(entry, chr1_reads_p, startpos, stoppos)        

                elif entry[0] == "2" and length >= min_s and length <= max_s:
                    load_read(entry, chr2_reads_p, startpos, stoppos)  

                elif entry[0] == "3" and length >= min_s and length <= max_s:    
                    load_read(entry, chr3_reads_p, startpos, stoppos)   

            elif entry[5] == "-":
                startpos = int(entry[2])
                stoppos = int(entry[1])
                if entry[0] == "1" and length >= min_s and length <= max_s:
                    load_read(entry, chr1_reads_m, startpos, stoppos)          

                elif entry[0] == "2" and length >= min_s and length <= max_s:
                    load_read(entry, chr2_reads_m, startpos, stoppos)

                elif entry[0] == "3" and length >= min_s and length <= max_s:   
                    load_read(entry, chr3_reads_m, startpos, stoppos)  
            else: 
                print("No matching strand specificity")
                # print(f'{startpos}:\t{chr3_reads.get(startpos)}')
                
    # save the objects as pickle objects to increase script runtime             
    save_pickles()

print(f'Collect reads between {min_s}nt and {max_s}nt size.') 
print(f'Collected:\nchr1:\t{len(chr1_reads_p.keys())} (+) and {len(chr1_reads_m)} (-)\nchr2:\t{len(chr2_reads_p)} (+) and {len(chr2_reads_m)} (-) \nchr3:\t{len(chr3_reads_p)} (+) and {len(chr3_reads_m)} (-)\n piRNAs.')            

# slide over all Flavi EVEs
print("Screen for piRNAs on EVEs...")

for EVE in EVEs_dict.keys():
    entry = EVEs_dict.get(EVE)
    chunk_num = 0
    pi_nums_p = 0                 # number of distinct piRNAs == individual 5' mapping positions
    pi_nums_m = 0                 # number of distinct piRNAs == individual 5' mapping positions
    # Structure of chunk_depth: {pos:[pos_mapping; neg_mapping]}
    chunk_depth_p = OrderedDict()           # piRNA quantification; how many are starting per position in chunk
    chunk_depth_m = OrderedDict()           # piRNA quantification; how many are starting per position in chunk

    if int(entry[3]) >= w_size:                               # if EVE longer than maximal subseq size for cloning
        ### positive orientation
        print(f'{EVE}: {entry[1]} - {int(entry[2])+1}')
        for num in range(int(entry[1]),int(entry[2])+1, 1):   # creating chunks, separated by minimal piRNA size
            chunk_num += 1
            if num+w_size <= int(entry[2]):                     # not over the right boundary!
                for pos in range(num, num+w_size):            # traverse the chunks
                    #chunk_lim = num+w_size                          
                    if int(entry[0]) == 1:                       
                        count_pi_perpos(chr1_reads_p, pos, num+w_size)
                        count_pi_perneg(chr1_reads_m, pos, num-w_size)
                    elif  int(entry[0]) == 2:
                        count_pi_perpos(chr2_reads_p, pos, num+w_size)
                        count_pi_perneg(chr2_reads_m, pos, num-w_size)
                    elif  int(entry[0]) == 3:
                        count_pi_perpos(chr3_reads_p, pos, num+w_size)
                        count_pi_perneg(chr3_reads_m, pos, num-w_size)
                chunk_name = f'chunk_+_{chunk_num}'
                
                chunkstart = num
                chunkstop = num+w_size-1
                # print(len(chunk_depth_p))
                # print(len(chunk_depth_m))
                sum_depth_p, mean_depth_p, max_depth_p = calc_metrics(pi_nums_p, chunk_depth_p)
                sum_depth_m, mean_depth_m, max_depth_m = calc_metrics(pi_nums_m, chunk_depth_m)

                
                # if pi_nums_p != 0:                    # if piRNAs found on this chunk, calculate mean and max depth
                #     depths_+ = list(chunk_depth_p.values())
                #     depths_+ = np.array(depths)
                #     mean_depth_p = np.mean(depths)
                #     max_depth_m = np.amax(depths)
                #     sum_depth_p = np.sum(depths)
                #     # print(f'EVE\tchunk_name\tchunkstart\tchunkstop\tmean_depth\tmax_depth\tpi_nums')                    
                # else:
                #     mean_depth_p = 0
                #     max_depth_m = 0   
                #     sum_depth_p = 0 
                
                # save the chunk to its EVE entry 
                # print(f'{EVE}\t{chunk_name}\t{chunkstart}\t{chunkstop}\t{mean_depth}\t{max_depth}\t{pi_nums}')
                if piRNA_EVE_p.get(EVE) != None:
                    piRNA_EVE_p[EVE].update({chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]})       # apend existing dict entry
                else:
                    piRNA_EVE_p.update({EVE:{chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]}})    # start new dict entry for new EVE
                
                pi_nums_p = 0
                chunk_depth_p = OrderedDict()
                pi_nums_m = 0
                chunk_depth_m = OrderedDict()

        # chunk_num = 0        

        ### negative orientation   
        # print(f'{EVE}: {entry[2]} - {int(entry[1])-1}')     
        # for num2 in range(int(entry[2]),int(entry[1])-1, -1):   # creating chunks, separated by minimal piRNA size
        #     chunk_num += 1
        #     # print(f'{num2} : {num2-w_size}')
        #     if num2-w_size >= int(entry[1]):                     # not over the right boundary!
        #         for neg in range(num2, num2-w_size, -1):            # traverse the chunks      
        #             # print(neg)                    
        #             if int(entry[0]) == 1:                       
        #                 count_pi_perneg(chr1_reads_m, neg, num2-w_size)
        #                 count_pi_perpos(chr1_reads_p, neg, num2+w_size)
        #             elif  int(entry[0]) == 2:
        #                 count_pi_perneg(chr2_reads_m, neg, num2-w_size)
        #                 count_pi_perpos(chr2_reads_p, neg, num2+w_size)
        #             elif  int(entry[0]) == 3:
        #                 count_pi_perneg(chr3_reads_m, neg, num2-w_size)
        #                 count_pi_perpos(chr3_reads_p, neg, num2+w_size)

        #             chunk_name = f'chunk_-_{chunk_num}'

        #         chunkstart = num2
        #         chunkstop = num2-w_size

        #         sum_depth_p, mean_depth_p, max_depth_p = calc_metrics(pi_nums_p, chunk_depth_p)
        #         sum_depth_m, mean_depth_m, max_depth_m = calc_metrics(pi_nums_m, chunk_depth_m)

        #         # print(f'{EVE}\t{chunk_name}\t{chunkstart}\t{chunkstop}\t{mean_depth}\t{max_depth}\t{pi_nums}')
        #         if piRNA_EVE_m.get(EVE) != None:
        #             piRNA_EVE_m[EVE].update({chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]})       # apend existing dict entry
        #         else:
        #             piRNA_EVE_m.update({EVE:{chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]}})    # start new dict entry for new EVE
                
        #         pi_nums_p = 0
        #         chunk_depth_p = OrderedDict()
        #         pi_nums_m = 0
        #         chunk_depth_m = OrderedDict()
            
    else: # all EVEs shorter 500nt do not need to be chunked, just statistics:
        ### pos orientation
        print(f'{EVE}: {entry[1]} - {int(entry[2])+1}')
        for pos in range(int(entry[1]),int(entry[2])+1):            # traverse the chunks
            chunk_lim = int(entry[2])                        
            if int(entry[0]) == 1:                       
                count_pi_perpos(chr1_reads_p, pos, int(entry[2]))
                count_pi_perneg(chr1_reads_m, pos, int(entry[1]))

            elif  int(entry[0]) == 2:
                count_pi_perpos(chr2_reads_p, pos, int(entry[2]))
                count_pi_perneg(chr2_reads_m, pos, int(entry[1]))

            elif  int(entry[0]) == 3:
                count_pi_perpos(chr3_reads_p, pos, int(entry[2]))
                count_pi_perneg(chr3_reads_m, pos, int(entry[1]))

        # calculate chunk data  
        chunk_name = f'chunk_+_1'
        chunkstart = int(entry[1])
        chunkstop = int(entry[2])

        # if EVE == "Fla4":
            # print(chunk_name)
            # print(chunk_depth_p)
            # print(chunk_depth_m)
        sum_depth_p, mean_depth_p, max_depth_p = calc_metrics(pi_nums_p, chunk_depth_p)
        sum_depth_m, mean_depth_m, max_depth_m = calc_metrics(pi_nums_m, chunk_depth_m)
        # if EVE == "Fla14":
        #     print(chunk_depth_p)
        #     print(chunk_depth_m)
        # # save the chunk to its EVE entry 
        if piRNA_EVE_p.get(EVE) != None:
            piRNA_EVE_p[EVE].update({chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]})       # apend existing dict entry
        else:
            piRNA_EVE_p.update({EVE:{chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]}})    # start new dict entry for new EVE
                
        pi_nums_p = 0
        chunk_depth_p = OrderedDict()
        pi_nums_m = 0
        chunk_depth_m = OrderedDict()

        ### negative orientation      
        # print(f'{EVE}: {entry[2]} - {int(entry[1])-1}') 
        # for neg in range(int(entry[2]),int(entry[1])-1, -1):   # creating chunks, separated by minimal piRNA size
                  
        #     if int(entry[0]) == 1:                       
        #         count_pi_perneg(chr1_reads_m, neg, int(entry[1]))
        #         count_pi_perpos(chr1_reads_p, neg, int(entry[2]))

        #     elif  int(entry[0]) == 2:
        #         count_pi_perneg(chr2_reads_m, neg, int(entry[1]))
        #         count_pi_perpos(chr2_reads_p, neg, int(entry[2]))
                
        #     elif  int(entry[0]) == 3:
        #         count_pi_perneg(chr3_reads_m, neg, int(entry[1]))
        #         count_pi_perpos(chr3_reads_p, neg, int(entry[2]))


        # # calculate chunk data      
        # chunk_name = f'chunk_-_1'
        # chunkstart = int(entry[2])
        # chunkstop = int(entry[1])
        # # if EVE == "Fla4":
        # #     print(chunk_name)
        # #     print(chunk_depth_p)
        # #     print(chunk_depth_m)
        # sum_depth_p, mean_depth_p, max_depth_p = calc_metrics(pi_nums_p, chunk_depth_p)
        # sum_depth_m, mean_depth_m, max_depth_m = calc_metrics(pi_nums_m, chunk_depth_m)
        # # print(f'{EVE}\t{chunk_name}\t{chunkstart}\t{chunkstop}\t{mean_depth}\t{max_depth}\t{pi_nums}')
        

        # if piRNA_EVE_m.get(EVE) != None:
        #     piRNA_EVE_m[EVE].update({chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]})       # apend existing dict entry
        # else:
        #     piRNA_EVE_m.update({EVE:{chunk_name:[chunkstart, chunkstop, pi_nums_p, max_depth_p, mean_depth_p, sum_depth_p,  chunk_depth_p, pi_nums_m, max_depth_m, mean_depth_m, sum_depth_m,  chunk_depth_m]}})    # start new dict entry for new EVE
                
        # pi_nums_p = 0
        # chunk_depth_p = OrderedDict()
        # pi_nums_m = 0
        # chunk_depth_m = OrderedDict()

# write chunk summary file:

for EVE in EVEs_dict.keys():
    with open(f'{outpath}/{EVE}_chunks_summary.csv', 'w') as chunkf:
        # print(f'{EVE}')
        data_pos = piRNA_EVE_p.get(EVE)
        data_neg = piRNA_EVE_m.get(EVE)
        
        chunkf.write("EVE SEQUENCE\n")
        chunkf.write(f'EVE_Name\tChromosome\tStart_Pos\tStop_Pos\tLength\tStrand\n')
        chunkf.write(f'>{EVE}\t{EVEs_dict[EVE][0]}\t{EVEs_dict[EVE][1]}\t{EVEs_dict[EVE][2]}\t{EVEs_dict[EVE][3]}\t{EVEs_dict[EVE][4]}\n')
        chunkf.write("CHUNK SEQUENCES\n")
        chunkf.write(f'Chunk_Name\tRead_Orientation\tStart_Pos\tStop_Pos\tpiRNA_Counts\tSum_Counts\tMean_Depth\tMax_Depth\n')
        # write pos orient
        for chunkp in data_pos.keys():
            # print(f'{chunkp}')
            chunk_data = data_pos.get(chunkp)
            # if EVE == "Fla2":
            #     print(chunk_data[4])
            #     print(chunk_data[9])
            chunkf.write(f'{chunkp}\t+\t{chunk_data[0]}\t{chunk_data[1]}\t{chunk_data[2]}\t{chunk_data[5]}\t{chunk_data[4]}\t{chunk_data[3]}\n')
            chunkf.write(f'{chunkp}\t-\t{chunk_data[0]}\t{chunk_data[1]}\t{chunk_data[7]}\t{chunk_data[10]}\t{chunk_data[9]}\t{chunk_data[8]}\n')
            # write piRNA counts per position of each chunk
            with open(f'{outpath}/raw_out/{EVE}_{chunkp}_raw.csv', 'w') as chunkraw_p:
                chunkraw_p.write(f'>{EVE}\t{chunkp}\t{chunk_data[0]}\t{chunk_data[1]}\n')
                chunkraw_p.write(f'Position\tpiRNA_counts_+\tpiRNA_counts_-\n')
                for pos in chunk_data[6].keys():
                    chunkraw_p.write(f'{pos}\t{chunk_data[6][pos]}\t{chunk_data[11][pos]}\n')
















        # write neg orient            
        # for chunkm in data_neg.keys():
        #     # print(f'{chunkm}')
        #     chunk_data = data_neg.get(chunkm)
        #     chunkf.write(f'{chunkm}\t+\t{chunk_data[0]}\t{chunk_data[1]}\t{chunk_data[2]}\t{chunk_data[5]}\t{chunk_data[4]}\t{chunk_data[3]}\n')
        #     chunkf.write(f'{chunkm}\t-\t{chunk_data[0]}\t{chunk_data[1]}\t{chunk_data[7]}\t{chunk_data[10]}\t{chunk_data[9]}\t{chunk_data[8]}\n')
        #     # write piRNA counts per position of each chunk
        #     with open(f'{outpath}/raw_out/{EVE}_{chunkm}_raw.csv', 'w') as chunkraw_m:
        #         chunkraw_m.write(f'>{EVE}\t{chunkm}\t{chunk_data[0]}\t{chunk_data[1]}\n')
        #         chunkraw_m.write(f'Position\tpiRNA_counts_+\tpiRNA_counts_-\n')
        #         for pos in chunk_data[6].keys():
        #             chunkraw_m.write(f'{pos}\t{chunk_data[6][pos]}\t{chunk_data[11][pos]}\n')
