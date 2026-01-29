# CiFi

Toolkit for downstream processing of CiFi long reads.

https://dennislab.org/cifi

## Install

```bash
pip install cifi
# or
mamba install bioconda::cifi
```

## Commands

| Command | Description |
|---------|-------------|
| `cifi qc` | Sample reads and report enzyme site frequency, fragment sizes, estimated yield |
| `cifi digest` | In-silico digestion → paired-end FASTQ (all pairwise contacts) |
| `cifi filter` | MAPQ-based filtering of aligned paired-end BAM |
| `cifi enzymes` | List built-in restriction enzymes |

### qc

```bash
cifi qc reads.bam -e HindIII -o qc_out
cifi qc reads.bam -e NlaIII -n 50000 -o qc_out    # sample 50k reads
cifi qc reads.bam -e HindIII -n 0 -o qc_out        # all reads
cifi qc reads.bam --site GANTC --cut-pos 1 -o qc_out  # custom site
```

Writes an output directory with HTML report, JSON, TSV tables, distribution plots (PNG), and a multi-page PDF.

### digest

```bash
cifi digest reads.bam -e HindIII -o output
cifi digest reads.fq.gz -e NlaIII -o output -m 5 --gzip
cifi digest reads.bam --site GANTC --cut-pos 1 -o output
```

Produces `{prefix}_R1.fastq` and `{prefix}_R2.fastq` (optionally gzipped), plus an HTML report and JSON stats.

### filter

```bash
cifi filter aligned.bam -o filtered.bam -q 30
```

Keeps properly paired reads where both mates meet the MAPQ threshold.

## Enzymes

Built-in enzymes:

| 4-cutters | 6-cutters |
|-----------|-----------|
| NlaIII (CATG) | HindIII (AAGCTT) |
| DpnII (GATC) | |
| MboI (GATC) | |
| Sau3AI (GATC) | |

Any recognition site can be specified with `--site` and `--cut-pos`, including IUPAC degenerate bases (N, R, Y, W, S, M, K, B, D, H, V).

## How it works

CiFi reads are concatemers of restriction fragments from genomic regions in 3D proximity. The toolkit finds all enzyme cut sites in each read, extracts fragments, and generates every pairwise combination as a pseudo paired-end read:

```
Read with 4 fragments: [A]-[B]-[C]-[D]
Pairs: A-B, A-C, A-D, B-C, B-D, C-D  (6 pairs)
```

## Citation

Coming soon.
