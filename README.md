# pyBlast

## Project Overview

**pyBlast** is a powerful Python package providing a command-line interface for BLAST (Basic Local Alignment Search Tool). BLAST is a fundamental tool in bioinformatics for comparing biological sequences, including DNA, RNA, and proteins. pyBlast aims to streamline the process of constructing and executing BLAST queries, offering an intuitive interface and extensive parameter options to support various complex analysis needs.

## Key Features

- **Extensive Parameter Options**: Supports a comprehensive range of BLAST parameters, including database selection, alignment tasks, and output formats.
- **Colored Help Information**: Utilizes ANSI color codes in help documentation for enhanced clarity in parameter explanations.
- **Local and Remote Execution**: Enables running BLAST searches both locally and remotely.
- **Multi-threading Support**: Accelerates searches using multiple threads for improved performance.
- **Filtering and Masking Options**: Includes various filtering and masking techniques such as DUST, WindowMasker, and custom strategies.
- **Support for Multiple BLAST Tools**: Includes `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, and `makeblastdb` for diverse sequence alignment needs.

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

`pyBlast` provides a set of command-line tools for performing BLAST searches. Use the `--help` flag with any tool to view supported options and parameters:

```bash
blastn --help
blastp --help
blastx --help
tblastn --help
tblastx --help
makeblastdb --help
```

### Tool Descriptions

- **`blastn`**: Aligns nucleotide sequences against a nucleotide database.
- **`blastp`**: Aligns protein sequences against a protein database.
- **`blastx`**: Translates nucleotide sequences into proteins and aligns them against a protein database.
- **`tblastn`**: Aligns protein sequences against a nucleotide database by translating the nucleotide sequences into proteins.
- **`tblastx`**: Translates both nucleotide sequences (query and database) into proteins and aligns them.
- **`makeblastdb`**: Converts raw sequence files into a searchable BLAST database.

### Example Commands

#### Perform a BLAST Search

```bash
blastn -query input_sequence.fasta -db nt -out results.txt -evalue 1e-5 -num_threads 4
```

#### Perform BLAST Search with a Specific Template Type

```bash
blastn -query input_sequence.fasta -db nt -out results.txt -template_type discontiguous -template_length 20
```

#### Create a BLAST Database

```bash
makeblastdb -in sequences.fasta -dbtype nucl -title "My Database" -out my_database
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

**pyBlast** is an open-source project dedicated to advancing bioinformatics tools. We encourage all developers and researchers interested in bioinformatics to contribute to this project.
