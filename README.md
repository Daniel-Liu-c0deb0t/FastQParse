# FastQParse
Java toolkit for preprocessing .fastq files for the GBS Pipeline

Download FastQParse.jar for the latest compiled binaries.

### Description:
FastQParse handles the initial processing of sequenced DNA for the GBS pipeline, which is necessary for base calling, alignment, and downstream analysis. FastQParse integrates many steps in needed when preprocessing into one single command. Single end, paired end, and UMI fastq files are all supported. Different matching algorithms are implemented to maximize the flexibility of the toolkit. Running different features separately allows FastQParse to handle files for other pipelines.

### Features:
- Demultiplex single or paired-end reads (two matching algorithms)
- Deduplicate (UMI)
- Support for UMI index files
- Merge paired-end reads (two matching algorithms)
- Quality trim (two trimming algorithms)
- Remove adapters (5' or 3' and anchored or not anchored, two matching algorithms)
- Gzip input/output
- Support for different types of mismatches (Hamming, Levenshtein, or Bayesian probability)
- Remove reads with low average quality or high expected error rate
- Comprehensive stats when demultiplexing (.stats file)
- Dump removed reads to other files
- Do all of the above in one command (can run demultiplexing, quality trimming, adapters trimming, merging paired ends, etc. together)
- Or use quality trim, remove adapters, merge paired ends, deduplicate, etc. as standalone features (ex. only remove adapters in one file)
- Parallel processing using multi-threading
- Efficient and flexible matching using bit-parallel algorithms and probability based matching
- And more...

Tutorial available [here](https://github.com/Daniel-Liu-c0deb0t/FastQParse/wiki/Tutorial).

List of commands available [here](https://github.com/Daniel-Liu-c0deb0t/FastQParse/wiki/Commands).

Algorithms and implementation details available [here](https://github.com/Daniel-Liu-c0deb0t/FastQParse/wiki/Algorithms).

Wiki and other information available [here](https://github.com/Daniel-Liu-c0deb0t/FastQParse/wiki).
