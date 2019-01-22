# Iterative-HMMER
Iterative HMMER search algorithm

# Requirements :
## Python libraries
* pandas
* numpy
* Bio
* matplotlib
* progressbar
## Software
* HMMER suite
## Protein Database
You will need a protein database to search against. We recommend you use TREMBL. The database needs to be in the `FASTA` format.

You can download TREMBL from this ftp: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

Informations here: https://www.uniprot.org/downloads
## Alignement
You will need an alignement to draw from such that we can build HMM with the rolling window. The alignment must be in `FASTA`.

# Usage

```Usage: Iterative_hmmer_search.py [options]

Options:
  -h, --help            show this help message and exit
  -a ALIGNPATH, --alignment=ALIGNPATH
                        [Required] Location of the FASTA alignement file.
  -w WINDOW, --window=WINDOW
                        [Required] Size of the sliding window (aka minimum nb
                        of columns to be used for HMM generation)
  -n NAME, --name=NAME  [Required] Name of the analysis)
  -l, --analyze         [Optional] Analyze only)
  -d DB, --db=DB        Path to DB (default to Trembl)
  -o, --oskar           Is it the Oskar Gene ?
  -p, --preprocess      Do not preprocess the search. (slower and more prone
                        to false positives)
```
