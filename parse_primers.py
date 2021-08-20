#!/home/max/bin/anaconda3/bin/python


import sys, re

infile = "/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/primer3/primers_260nt.csv"
outfile = "/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/primer3/primers_order.csv"
primer_fwd = []
primer_rev = []

count_l = 0
count_r = 0


with open(infile , 'r') as csv:
    for line in csv:
        if line.startswith("SI"):
            entry = line.split("\t")
            print(entry)
            target = entry[0]
            primer = entry[1]
            seq = entry[10]
            print(target)
            print()
            if bool(re.match(".*LEFT.*", entry[1])):

                count_l += 1
                primer_fwd.append([target, count_l, primer, seq])
            elif bool(re.match(".*RIGHT.*", entry[1])):
                count_r += 1
                primer_rev.append([target, count_r, primer, seq])

with open(outfile, 'w') as out:
    for record in primer_fwd:
        out.write(f'SINV_GFP_fwd_{record[1]}\t{record[3]}\t{record[0]}\t{record[2]}\n')
    for record in primer_rev:
        print(record)
        out.write(f'SINV_GFP_rev_{record[1]}\t{record[3]}\t{record[0]}\t{record[2]}\n')