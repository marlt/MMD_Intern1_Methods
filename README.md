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