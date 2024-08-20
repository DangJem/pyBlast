# pyBlast

## Project Overview

**pyBlast** is a versatile Python package that provides a command-line interface for BLAST (Basic Local Alignment Search Tool). BLAST is a widely used tool in bioinformatics for comparing biological sequences, including DNA, RNA, and proteins. `pyBlast` offers an efficient way to construct and execute BLAST queries with a rich set of parameters to accommodate diverse analysis needs.

## Key Features

- **Extensive Parameter Options**: Supports a broad range of BLAST parameters including database selection, alignment tasks, and output formats.
- **Colored Help Information**: Uses ANSI color codes in help documentation for clearer parameter explanations.
- **Local Execution**: Designed to be run locally from downloaded binaries.
- **Multi-threading Support**: Utilizes multiple threads to enhance performance.
- **Filtering and Masking Options**: Offers various filtering and masking strategies like DUST and WindowMasker.
- **Support for Multiple BLAST Tools**: Includes `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `makeblastdb`, and `blast_formatter`.

## Installation

### Download and Install

1. **Download the Latest Release**: Visit the [releases page](https://github.com/DangJem/pyblast/releases) to download the latest release of `pyBlast`, for example, `pyBlast-0.0.6.zip`.

2. **Extract the Downloaded Archive**:
   
   ```bash
   unzip pyBlast-0.0.6.zip
   ```

3. **Navigate to the Extracted Directory**:
   
   ```bash
   cd pyBlast-0.0.6
   ```

4. **Run the Software**: Directly use the software from this directory or move it to a location included in your `PATH`.

## Usage Instructions

### Basic Usage

`pyBlast` provides several command-line tools for performing BLAST searches. Use the `-h` or `--help` flag with any tool to view supported options and parameters:

```bash
python pyBlast.py blastn -h
python pyBlast.py blastp -h
python pyBlast.py blastx -h
python pyBlast.py tblastn -h
python pyBlast.py tblastx -h
python pyBlast.py makeblastdb -h
python pyBlast.py blast_formatter -h
```

### Tool Descriptions

- **`blastn`**: Performs nucleotide-nucleotide sequence alignments. Aligns nucleotide sequences against a nucleotide database to find similarities. Useful for identifying homologous sequences, discovering genes, and annotating genome data.
    - Key Features:
        - High-speed searches for nucleotide sequences against nucleotide databases.
        - Various filtering and scoring options to refine search results.
        - Output formats for easy integration with other tools and analysis pipelines.
    - Typical Usage:
        - `python pyBlast.py blastn -query query.fasta -db database -out results.txt -evalue 1e-10`
        - This command performs a BLAST search using sequences from `query.fasta` against the database `database`, with results saved to `results.txt` and an E-value threshold of `1e-10`.
    - Additional Options:
        - `-query`: Input file containing nucleotide sequences to search.
        - `-db`: Name of the nucleotide database to search against.
        - `-out`: File to write the results to.
        - `-evalue`: E-value threshold for reporting matches.
        - `-task`: Specific BLASTN task to use (e.g., megablast, dc-megablast).

- **`makeblastdb`**: Creates a BLAST database from input sequence files. Converts raw sequence data into a searchable database format for BLAST. Supports nucleotide and protein sequences.
    - Key Features:
        - Converts raw sequence data into a searchable database format.
        - Supports nucleotide and protein sequences.
        - Allows for custom database titles and descriptions.
    - Typical Usage:
        - `python pyBlast.py makeblastdb -in sequences.fasta -dbtype nucl -out my_database`
        - This command reads sequences from `sequences.fasta` and creates a nucleotide database named `my_database`.
    - Additional Options:
        - `-in`: Input file containing the sequences.
        - `-dbtype`: Type of database (nucl for nucleotide, prot for protein).
        - `-out`: Output database file name.
        - `-title`: Custom title for the database.
        - `-taxid`: Taxonomy ID for the sequences, if applicable.

- **`blastp`**: Performs protein-protein sequence alignments. Compares an input protein sequence against a protein database to find homologous proteins and infer their functional and evolutionary relationships.
    - Key Features:
        - Searches for protein matches in a protein database.
        - Provides detailed alignment results including percent identity and alignment scores.
        - Flexible output formats for integration into various analysis workflows.
    - Typical Usage:
        - `python pyBlast.py blastp -query protein_query.fasta -db protein_db -out results.txt -evalue 1e-5`
        - This command aligns sequences from `protein_query.fasta` against `protein_db`, saving results to `results.txt` with an E-value threshold of `1e-5`.
    - Additional Options:
        - `-query`: File containing protein sequences to align.
        - `-db`: Name of the protein database to search against.
        - `-out`: File to write the results to.
        - `-evalue`: E-value threshold for reporting matches.
        - `-matrix`: Scoring matrix to use for alignment (e.g., BLOSUM62).

- **`blastx`**: Performs protein query vs. translated nucleotide database search. Aligns protein queries against a nucleotide database that is dynamically translated in all reading frames, allowing for comparisons across species and genomic analysis.
    - Key Features:
        - Translate nucleotide sequences dynamically for protein-level comparisons.
        - Comprehensive filtering and scoring mechanisms to refine your results.
        - Customizable output formats for various downstream applications and analyses.
    - Typical Usage:
        - `python pyBlast.py blastx -query query.fasta -db protein_db -out results.txt -evalue 1e-5`
        - This command performs a BLASTX search using sequences from `query.fasta` against the protein database `protein_db`, with results saved to `results.txt` and an E-value threshold of `1e-5`.
    - Additional Options:
        - `-query`: Input file containing nucleotide sequences to be translated and searched.
        - `-db`: Name of the protein database to search against.
        - `-out`: File to write the results to.
        - `-evalue`: E-value threshold for reporting matches.
        - `-task`: Specific BLASTX task to use (e.g., blastx-fast, blastx).

- **`tblastn`**: Performs protein query against translated nucleotide database search. Aligns protein sequences against a nucleotide database that is dynamically translated into all reading frames.
    - Key Features:
        - Translates nucleotide databases into protein sequences dynamically for alignment against protein queries.
        - Various filtering, scoring, and alignment customization options to refine search results.
        - Supports different output formats for downstream analysis and integration with other bioinformatics tools.
    - Typical Usage:
        - `python pyBlast.py tblastn -query protein.fasta -db nucleotide_db -out results.txt -evalue 1e-3`
        - This command performs a TBLASTN search using sequences from `protein.fasta` against the nucleotide database `nucleotide_db`, with results saved to `results.txt` and an E-value threshold of `1e-3`.
    - Additional Options:
        - `-query`: Input file containing protein sequences to search against the translated nucleotide database.
        - `-db`: Name of the nucleotide database to search, which will be dynamically translated into protein sequences.
        - `-out`: File to write the results to.
        - `-evalue`: E-value threshold for reporting significant matches.
        - `-task`: Specific TBLASTN task to execute (e.g., tblastn, tblastn-fast).

- **`tblastx`**: Performs translated nucleotide-nucleotide sequence alignments. Aligns nucleotide sequences against a nucleotide database by translating both query and database sequences into proteins.
    - Key Features:
        - Translated nucleotide-to-nucleotide alignments.
        - Various filtering and scoring options to enhance search results.
        - Output formats for integration with other tools and pipelines.
    - Typical Usage:
        - `python pyBlast.py tblastx -query query.fasta -db database -out results.txt -evalue 1e-10`
        - This command performs a BLAST search using translated sequences from `query.fasta` against the translated nucleotide database `database`, with results saved to `results.txt` and an E-value threshold of `1e-10`.
    - Additional Options:
        - `-db`: Name of the nucleotide database to search against.
        - `-query`: Input file containing nucleotide sequences to translate and search.
        - `-out`: File to write the search results.
        - `-evalue`: E-value threshold for reporting matches.
        - `-word_size`: Word size for the wordfinder algorithm.
        - `-qcov_hsp_perc`: Percent query coverage per HSP.
        - `-max_hsps`: Maximum number of HSPs per subject sequence to save for each query.
        - `-xdrop_ungap`: X-dropoff value for ungapped extensions.
        - `-searchsp`: Effective length of the search space.
        - `-sum_stats`: Enable or disable sum statistics.
        - `-max_intron_length`: Maximum intron length.
        - `-seg`: SEG filtering options.
        - `-soft_masking`: Apply filtering locations as soft masks.
        - `-matrix`: Scoring matrix name.
        - `-threshold`: Minimum score threshold for extending hits.
        - `-culling_limit`: Limit on hits enveloped by higher-scoring hits.
        - `-best_hit_overhang`: Best hit algorithm overhang value.
        - `-best_hit_score_edge`: Best hit algorithm score edge value.
        - `-subject_besthit`: Enable the best hit per subject sequence.
        - `-window_size`: Window size for multiple hits (use 0 for 1-hit algorithm).
        - `-lcase_masking`: Enable lowercase masking in query and subject sequences.
        - `-query_loc`: Location on the query sequence in 1-based offsets (Format: start-stop).
        - `-strand`: Strand to search against (e.g., plus, minus).
        - `-parse_deflines`: Enable parsing of query and subject deflines.
        - `-query_gencode`: Genetic code to use for query sequence translation.
        - `-db_gencode`: Genetic code to use for database sequence translation.
        - `-outfmt`: Output alignment format.
        - `-show_gis`: Display NCBI GI numbers in deflines.
        - `-num_descriptions`: Number of descriptions to show in the output.
        - `-num_alignments`: Number of alignments to show in the output.
        - `-line_length`: Line length for the output.
        - `-html`: Generate output in HTML format.
        - `-sorthits`: How to sort hits in the output.
        - `-sorthsps`: How to sort HSPs in the output.
        - `-max_target_seqs`: Maximum number of target sequences to report.
        - `-num_threads`: Number of threads to use for parallel processing.
        - `-remote`: Perform the search using the NCBI BLAST service instead of local resources.


### Example Commands

#### Perform a BLAST Search

```bash
python pyBlast.py blastn -query my_query.fasta -db my_database -out results.out
```

#### Create a BLAST Database

```bash
python pyBlast.py makeblastdb -in my_sequences.fasta -dbtype nucl -title "My Nucleotide Database" -out my_database
```

#### Perform Protein Alignments

```bash
python pyBlast.py blastp -query my_protein_query.fasta -db my_protein_database -out results.out
```

#### Perform Protein Alignments with Translated Nucleotide Sequences

```bash
python pyBlast.py blastx -query my_nucleotide_query.fasta -db my_protein_database -out results.out
```

#### Perform Protein-Nucleotide Alignments

```bash
python pyBlast.py tblastn -query my_protein_query.fasta -db my_nucleotide_database -out results.out
```

#### Perform Translated Nucleotide Alignments

```bash
python pyBlast.py tblastx -query my_nucleotide_query.fasta -db my_nucleotide_database -out results.out
```

For detailed help on each subcommand, use:

```bash
python pyBlast.py blastn -h
python pyBlast.py makeblastdb -h
python pyBlast.py blastp -h
python pyBlast.py blastx -h
python pyBlast.py tblastn -h
python pyBlast.py tblastx -h
```

## Contributing

We welcome contributions! If you have bug reports, feature requests, or pull requests, please follow these steps:

1. **Fork the Repository**: Click the "Fork" button on the repository page to create a copy under your own GitHub account.
2. **Clone Your Fork**: Clone the repository to your local machine.
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
    git add.
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
