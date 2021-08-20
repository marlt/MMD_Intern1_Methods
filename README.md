# Online Methods

This is the online methods repository for the MMD Internship 1 of Maximilian Arlt to the topic:

## SINV Primer Design
The primers were designed using the primer3 (v2.5.0) CLI on a Debian Buster operating system. It was performed in 5 rounds (+ last chunk), based on an individual settings file. The parameters were losened every round arbitrarily, keeping them "as small as possible".

The output primer choices got manually assessed on the best parameters to choose one pair per chunk. Chunks were the Sindbis virus genome chopped in 260 nt (with 10nt overlap each)

The programm is called via 

```
primer3 --strict_tags --p3_settings_file=/keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/primer3/settings/primer3_settings_round2.txt /keller2/Studium/Masterstudium/MMD/#studies/1year/Internship_1/data/primer3/input/GFP_SINV_chunks260_round2.fasta | ../scripts/count_hits.py
```

Output was transformed to make it easier to read:

```
../scripts/primers_csv.py ../data/primer3/output/GFP_SINV_chunks260_round2_out.txt ../data/primer3/output/GFP_SINV_chunks260_round2_out.csv
```

Final parameters were parsed into the right format using:

```
parser_csv.py 
```

## Flavi-nrEVE Chunk Selection
The script "subeveid.py" collects all read starting positions and the total counts separately for forward and reverse mapping reads

It screens all 500 nt chunks per nrEVE and compiles piRNA density and read counts as well as average reads per piRNA.
For this, it
    - screens fwd and rev mapping reads separately

    - counts only piRNAs that do not step over current window borders

    - screens the whole EVE, if it is shorter than 500 nt
    
    - was run considering every sRNA read between between 23 and 31 nt as piRNA

- output for every chunk: csv-file of the read depth per position in each chunk, separated for fwd and rev reads

- output for every nrEVE: compiled read sum, piRNA density and average coverage per piRNA for fwd and rev reads


script subeveid.py was run the following way

```

./subeveid.py 23 31 <outdir>

```
