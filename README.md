# CATHI
Implementation of the CATHI algorithm for identifying SiRTAs.

### Algorithm
This approach scores genomic sequences by counting strings of 4 or more nucleotides containing only consecutive `G` or `T`. While the `T`s must be single, there can be multiple `G`s in a row. However, strings may not consist of only `G`s.

By default, penalties are imposed for each occurrence of `GGTGG` or `TT` on the edges of a match. The default penalty is 1, but these can be changed by using the `-p`/`--penalty` option for `GGTGG` and the `-t`/`--ttpenalty` option for `TT`.

Sequences should be provided in FASTA format, where the header is minimally formatted with the chromosome and genomic coordinates (e.g. `CHR14:778653-778953 (-)`). The chromosome and coordinates can be separated by any of the following characters: `:`, `|`, or `-`.

There are two modes:
1. **Calculate the maximum score for the sequence across all windows.**
    - This is the default behavior.
    - For each sequence in the FASTA file, windows are generated (customized using the window (`-w`) and step size (`-s`) options) and the score is calculated for each window. Only the maximal score for each sequence is returned.
    - The script outputs the header of each sequence followed by the score.

Sample usage:
```
$ python cathi_score.py -w 100 -s 50 -p 1 -t 0 chr14.fa
```
This will calculate the max score for each sequence using a window size of 100bp, a step size of 50bp, a `GGTGG` penalty of 1, and a `TT` penalty of 0.

2. **Calculate and return a score for each window within each sequence.**
    - This mode is accessed using the `--signal` option.
    - For each sequence in the FASTA file, windows are generated (customized using the window (`-w`) and step size (`-s`) options) and the score is calculated for each window. The score for each window is returned.
    - A lower bound can be placed on the score by using the `--thresh` option. Only scores with a value >= the threshold will be returned. By default, the threshold is 0.
    - The strand of the sequence can be provided using the `--strand` option with either `+` or `-`. This will determine the value of the genomic coordinates returned. The default is `+`.
    - The output is provided in a four-column `bedgraph` format. The columns are `[chrom] [start] [end] [score]`, where the genomic coordinates represent the beginning and end of the window. The first line is the header of the sequence (preceded by a `#`).

Sample usage:
```
$ python cathi_score.py -w 100 --signal --thresh 20 chr14.fa
```
This will calculate the score for each sequence and window using a window size of 100bp, default step size, default penalties, and a threshold of 20. All windows returned will have a score of at least 20. If multiple sequences occur in the input file, output will be separated by a header line (e.g. `CHR14:778653-778953(-)`).


### Dependencies
This script requres [Biopython](https://biopython.org) and [NumPy](https://numpy.org). These can be easily installed using Anaconda.

```
conda install numpy
conda install -c conda-forge biopython
```

### Usage
```
usage: cathi_score.py [-h] [-p PENALTY] [-t TTPENALTY] [-w WINDOW] [-s STEP] [--signal] [--thresh THRESH] sequence_file

Calculate SiRTA score for sequence.

positional arguments:
  sequence_file         FASTA file of putative SiRTA sequences

optional arguments:
  -h, --help            show this help message and exit
  -p PENALTY, --penalty PENALTY
                        penalty to apply for each GGTGG occurence; default=1
  -t TTPENALTY, --ttpenalty TTPENALTY
                        penalty to apply for each flanking TT occurence; default=1
  -w WINDOW, --window WINDOW
                        sliding window size; default=100bp
  -s STEP, --step STEP  step size for sliding windows; default=1bp
  --signal              flag to print signal output instead of max score
  --thresh THRESH       when used with --signal, only return scores above this value; default=0
```
