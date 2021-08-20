#!/home/max/bin/anaconda3/bin/python

import sys, re
import collections

infile = sys.argv[1]    # primer3 output file
outfile = sys.argv[2]   # csv primer file
outfasta = '/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/primer3/input/GFP_SINV_chunks260_round7.fasta'
seqs_out = {}

headline=f'SEQ_NAME\tPRIMER\tSTART\tLEN\tTM\tGC\tSelfAnyTh\tSelfEndTh\tHpTh\tEND_STAB\tSEQ\tPROD_SIZE\tCompAnyTh\tCompEndTh\n'

pair = ["LEFT","RIGHT"]
nohit = 0
with open(outfasta, 'w') as out_fa:
    with open(outfile, 'w') as out_csv:
        out_csv.write(f'{headline}\n')
        with open(infile, 'r') as csv:
            entryline = ""
            record = {}
            for line in csv:
                entry = line.rstrip("\n").split("=")
                if not line.startswith("="):
                    record.update({entry[0]:entry[1]})
                else:
                    pairnum = int(record.get("PRIMER_PAIR_NUM_RETURNED"))
                    seq_name = record.get("SEQUENCE_ID")
                    print(seq_name)
                    if pairnum != 0:
                        for num in range(pairnum):
                            for side in pair:
                                loc = record.get(f'PRIMER_{side}_{num}')
                                loc = loc.split(",")
                                primname = f'PRIM_{side}_{num}'
                                start = loc[0]
                                length = loc[1]
                                TM = record.get(f'PRIMER_{side}_{num}_TM')
                                GC = record.get(f'PRIMER_{side}_{num}_GC_PERCENT')
                                SelfAnyTh = record.get(f'PRIMER_{side}_{num}_SELF_ANY_TH')
                                SelfEndTh = record.get(f'PRIMER_{side}_{num}_SELF_END_TH')
                                HpTh = record.get(f'PRIMER_{side}_{num}_HAIRPIN_TH')
                                E_Stab = record.get(f'PRIMER_{side}_{num}_END_STABILITY')
                                prim_seq = record.get(f'PRIMER_{side}_{num}_SEQUENCE')
                                size = record.get(f'PRIMER_PAIR_{num}_PRODUCT_SIZE')
                                comp_any_Th = record.get(f'PRIMER_PAIR_{num}_COMPL_ANY_TH')
                                comp_end_Th = record.get(f'PRIMER_PAIR_{num}_COMPL_END_TH')
                                out_csv.write(f'{seq_name}\t{primname}\t{start}\t{length}\t{TM}\t{GC}\t{SelfAnyTh}\t{SelfEndTh}\t{HpTh}\t{E_Stab}\t{prim_seq}\t{size}\t{comp_any_Th}\t{comp_end_Th}\n')
                            out_csv.write("\n")
                    else:
                        nohit += 1
                        noentry = record.get("PRIMER_PAIR_NUM_RETURNED")
                        out_fa.write(f'SEQUENCE_ID={seq_name}\n')
                        out_fa.write(f'SEQUENCE_TEMPLATE={record.get("SEQUENCE_TEMPLATE")}\n')
                        out_fa.write(f'=\n')
                        print(noentry)
print(f'Extracted seqs for next round: {nohit}')
                        # add to new primer entry file
                        # out.write(f'{entryline.rstrip("\t")}\n')


