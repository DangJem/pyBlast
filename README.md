# pyBlast

## Project Overview

**pyBlast** is a powerful Python package that provides a command-line interface for BLAST (Basic Local Alignment Search Tool). BLAST is a standard tool in bioinformatics used for comparing biological sequences. It can be used to compare DNA, RNA, and protein sequences. pyBlast offers a convenient way to construct and execute BLAST queries, supporting a wide range of parameters to accommodate complex analysis needs.

The goal of this project is to make performing BLAST searches in Python more efficient and provide a rich interface for customizing search parameters.

## Key Features

- **Extensive Parameter Options**: Supports nearly all BLAST parameters including database selection, alignment tasks, and output formats.
- **Colored Help Information**: Uses ANSI color codes in help documentation for clearer parameter explanations.
- **Remote and Local Execution**: Supports running BLAST locally and remotely via network execution.
- **Multi-threading Support**: Accelerates BLAST searches using multiple threads for improved performance.
- **Filtering and Masking Options**: Provides various filtering and masking options, such as DUST, WindowMasker, and custom masking strategies.

## Installation

### Install from PyPI

You can install `pyBlast` from PyPI using pip:

```bash
pip install pyblast
```

### Install from Source

1. Clone the project repository:

    ```bash
    git clone https://github.com/DangJem/pyblast.git
    ```

2. Navigate into the project directory:

    ```bash
    cd pyblast
    ```

3. Install the project:

    ```bash
    pip install .
    ```

## Usage Instructions

### Basic Usage

pyBlast provides a command-line tool `blastn` for executing BLAST searches. Use the `--help` flag to view all supported options and parameters:

```bash
blastn --help
```

### Parameter Description

Here are descriptions for the main parameters of `blastn`:

- **`-help`**: Print USAGE, DESCRIPTION, and ARGUMENTS; ignore all other parameters.
- **`-version`**: Print version number; ignore other arguments.
- **`-query`**: Input file name.
- **`-query_loc`**: Location on the query sequence in 1-based offsets (Format: start-stop).
- **`-strand`**: Query strand(s) to search against database/subject.
- **`-task`**: Task to execute.
- **`-db`**: BLAST database name.
- **`-out`**: Output file name.
- **`-evalue`**: Expectation value (E) threshold for saving hits.
- **`-word_size`**: Word size for the wordfinder algorithm (length of best perfect match).
- **`-gapopen`**: Cost to open a gap.
- **`-gapextend`**: Cost to extend a gap.
- **`-penalty`**: Penalty for a nucleotide mismatch.
- **`-reward`**: Reward for a nucleotide match.
- **`-use_index`**: Use MegaBLAST database index.
- **`-index_name`**: MegaBLAST database index name (deprecated; use only for old style indices).
- **`-subject`**: Subject sequence(s) to search.
- **`-subject_loc`**: Location on the subject sequence in 1-based offsets (Format: start-stop).
- **`-outfmt`**: Alignment view options.
- **`-show_gis`**: Show NCBI GIs in deflines?
- **`-num_descriptions`**: Number of database sequences to show one-line descriptions for.
- **`-num_alignments`**: Number of database sequences to show alignments for.
- **`-line_length`**: Line length for formatting alignments.
- **`-html`**: Produce HTML output?
- **`-sorthits`**: Sorting option for hits.
- **`-sorthsps`**: Sorting option for HSPs.
- **`-dust`**: Filter query sequence with DUST (Format: "yes", "level window linker", or "no" to disable).
- **`-filtering_db`**: BLAST database containing filtering elements (i.e., repeats).
- **`-window_masker_taxid`**: Enable WindowMasker filtering using a Taxonomic ID.
- **`-window_masker_db`**: Enable WindowMasker filtering using this repeats database.
- **`-soft_masking`**: Apply filtering locations as soft masks.
- **`-lcase_masking`**: Use lower case filtering in query and subject sequence(s)?
- **`-gilist`**: Restrict search of database to list of GIs.
- **`-seqidlist`**: Restrict search of database to list of SeqIDs.
- **`-negative_gilist`**: Restrict search of database to everything except the specified GIs.
- **`-negative_seqidlist`**: Restrict search of database to everything except the specified SeqIDs.
- **`-taxids`**: Restrict search of database to include only the specified taxonomy IDs and their descendants (multiple IDs delimited by ",").
- **`-negative_taxids`**: Restrict search of database to everything except the specified taxonomy IDs and their descendants (multiple IDs delimited by ",").
- **`-taxidlist`**: Restrict search of database to include only the specified taxonomy IDs and their descendants.
- **`-negative_taxidlist`**: Restrict search of database to everything except the specified taxonomy IDs and their descendants.
- **`-no_taxid_expansion`**: Do not expand the taxonomy IDs provided to their descendant taxonomy IDs.
- **`-entrez_query`**: Restrict search with the given Entrez query.
- **`-db_soft_mask`**: Filtering algorithm ID to apply to the BLAST database as soft masking.
- **`-db_hard_mask`**: Filtering algorithm ID to apply to the BLAST database as hard masking.
- **`-perc_identity`**: Percent identity.
- **`-qcov_hsp_perc`**: Percent query coverage per HSP.
- **`-max_hsps`**: Set maximum number of HSPs per subject sequence to save for each query.
- **`-culling_limit`**: If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit.
- **`-best_hit_overhang`**: Best Hit algorithm overhang value (recommended value: 0.1).
- **`-best_hit_score_edge`**: Best Hit algorithm score edge value (recommended value: 0.1).
- **`-subject_besthit`**: Turn on best hit per subject sequence.
- **`-max_target_seqs`**: Maximum number of aligned sequences to keep.
- **`-template_type`**: Discontiguous MegaBLAST template type.
- **`-template_length`**: Discontiguous MegaBLAST template length.
- **`-dbsize`**: Effective length of the database.
- **`-searchsp`**: Effective length of the search space.
- **`-import_search_strategy`**: Search strategy to use.
- **`-export_search_strategy`**: File name to record the search strategy used.
- **`-xdrop_ungap`**: X-dropoff value (in bits) for ungapped extensions.
- **`-xdrop_gap`**: X-dropoff value (in bits) for preliminary gapped extensions.
- **`-xdrop_gap_final`**: X-dropoff value (in bits) for final gapped alignment.
- **`-no_greedy`**: Use non-greedy dynamic programming extension.
- **`-min_raw_gapped_score`**: Minimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages.
- **`-ungapped`**: Perform ungapped alignment only?
- **`-window_size`**: Multiple hits window size, use 0 to specify 1-hit algorithm.
- **`-off_diagonal_range`**: Number of off-diagonals to search for the 2nd hit, use 0 to turn off.
- **`-parse_deflines`**: Should the query and subject defline(s) be parsed?
- **`-num_threads`**: Number of threads (CPUs) to use in the BLAST search.
- **`-mt_mode`**: Multi-thread mode to use in BLAST search.
- **`-remote`**: Execute search remotely?

## Examples

Here are some usage examples for `pyBlast`:

### Perform a BLAST Search

```bash
blastn -query input_sequence.fasta -db nt -out results.txt -evalue 1e-5 -num_threads 4
```

### Perform BLAST Search with a Specific Template Type

```bash
blastn -query input_sequence.fasta -db nt -out results.txt -template_type discontiguous -template_length 20
```

## Contributing

We welcome contributions! If you have any bug reports, feature requests, or pull requests, please follow these steps:

1. **Fork the Repository**: Click the "Fork" button on the repository page to create a copy of the repository under your own GitHub account.
2. **Clone

 Your Fork**: Clone the repository to your local machine.
    ```bash
    git clone https://github.com/yourusername/pyblast.git
    ```
3. **Create a Branch**: Create a new branch for your changes.
    ```bash
    git checkout -b feature/your-feature
    ```
4. **Make Changes**: Implement your changes.
5. **Commit Changes**: Commit your changes to your branch.
    ```bash
    git add .
    git commit -m "Describe your changes"
    ```
6. **Push Changes**: Push your branch to your fork.
    ```bash
    git push origin feature/your-feature
    ```
7. **Create a Pull Request**: Open a pull request on GitHub to propose your changes.

## License

This project is licensed under the [GNU General Public License (GPL)](LICENSE). See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or suggestions, please reach out:

- **Email**: dangjem0730@gmail.com
- **GitHub Issues**: [Open an Issue](https://github.com/DangJem/pyblast/issues)

Thank you for your interest and support!

---

**pyBlast** is an open-source project dedicated to advancing the use and development of bioinformatics tools. We encourage all developers and researchers interested in bioinformatics to contribute to this project.


### Changes Made:
1. **License Section**: Updated to reflect GPL licensing.
2. **General Cleanup**: Ensured the text is clear and consistent with the GPL license requirements.

Feel free to adjust any details or add additional sections as needed!
