# CATHI
Implementation of the CATHI algorithm for identifying SiRTAs.

### Algorithm
This approach scores genomic sequences by counting strings of 4 or more nucleotides containing only consecutive `G` or `T`. The sequence must always start with a `G`, and, while the `T`s must be single, there can be multiple `G`s in a row (up to 3). Strings may not consist of only `G`s, only `GGTGG` (and expansions, e.g. `GGTGGTGG`), or only `GTGGTGG`.

Penalties may be imposed for each occurrence of `GGTGG` or for a `TT` on the edge of a match (front or back). There are no penalties applied by default, but this can be changed by using the `-p`/`--penalty` option for `GGTGG` and the `-t`/`--ttpenalty` option for `TT`.

Sequences should be provided in FASTA format, where the header is minimally formatted with the chromosome and genomic coordinates (e.g. `CHR14:778653-778953 (-)`). The chromosome and coordinates can be separated by any of the following characters: `:`, `|`, or `-`.

There are two modes:
1. **Calculate the maximum score for the sequence across all windows.**
    - This is the default behavior.
    - For each sequence in the FASTA file, windows are generated (customized using the window (`-w`) and step size (`-s`) options) and the score is calculated for each window. Only the maximal score for each sequence is returned.
    - By default, the script outputs the header of each sequence followed by the score. If you prefer BED formatted output, use the `--bedformat` option to print the location followed by the maximal score for each sequence. (This option assumes the header contains minimally formatted with the chromosome name and start coordinate: e.g., `CHR14:778653` or `14-268677`).

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
    - If you would like to merge any overlapping coordinates in the output `bedgraph`, use the `--cluster` option. This will merge all overlapping (and directly adjacent) windows into a single region where the beginning is the start coordinate of the most upstream window and the ending is end coordinate of the most downstream window. The score is the maximum CATHI value across all merged windows.

Sample usage:
```
$ python cathi_score.py -w 100 --signal --thresh 20 chr14.fa
```
This will calculate the score for each sequence and window using a window size of 100bp, default step size, default penalties, and a threshold of 20. All windows returned will have a score of at least 20. If multiple sequences occur in the input file, output will be separated by a header line (e.g. `CHR14:778653-778953(-)`).


### Dependencies
This script requres several common biocomputing packages: [Biopython](https://biopython.org),  [NumPy](https://numpy.org), [Pandas](https://pandas.pydata.org), and [pybedtools](https://daler.github.io/pybedtools/index.html#). These can be easily installed using Anaconda. 

```
conda install numpy pandas regex
conda install -c conda-forge biopython
conda install -c conda-forge -c bioconda pybedtools
```
Please see the links for more detailed installation instructions for each package (including `pip` and installation from source).

### Usage
```
usage: cathi_score.py [-h] [-p PENALTY] [-t TTPENALTY] [-w WINDOW] [-s STEP] [--bedformat] [--signal]
                      [--thresh THRESH] [--strand {+,-}] [--cluster] sequence_file

Calculate SiRTA score for sequence.

positional arguments:
  sequence_file         FASTA file of putative SiRTA sequences

optional arguments:
  -h, --help            show this help message and exit
  -p PENALTY, --penalty PENALTY
                        penalty to apply for each GGTGG occurence; default=0
  -t TTPENALTY, --ttpenalty TTPENALTY
                        penalty to apply for each flanking TT occurence; default=0
  -w WINDOW, --window WINDOW
                        sliding window size; default=100bp
  -s STEP, --step STEP  step size for sliding windows; default=1bp
  --bedformat           print max score for seq (across windows) with location in BED format; default=False
  --signal              flag to print signal output instead of max score
  --thresh THRESH       when used with --signal, only return scores above this value; default=0
  --strand {+,-}        flag to specify strand with --signal; default=+
  --cluster             flag to merge overlapping windows in output in --signal mode; default=Fals
```
