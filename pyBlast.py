import argparse
import subprocess
import shutil
import sys

class BlastTool:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
        description='BLAST Tools for sequence alignment. This suite includes tools for sequence alignment with BLAST, '
                    'database creation, and querying functionalities.',
        formatter_class=argparse.RawTextHelpFormatter)
        self.subparsers = self.parser.add_subparsers(dest='command')


self.parser.epilog = (
    '\033[1;34mBLAST (Basic Local Alignment Search Tool)\033[0m is a widely used bioinformatics tool designed to identify similarities between sequences. '
    'It helps in functional, structural, and evolutionary studies of sequences by finding local alignments between them.\n\n'
    'The suite includes:\n'
    '  \033[1;32m1. blastn:\033[0m Perform nucleotide-nucleotide alignments\n'
    '     Aligns nucleotide sequences against a nucleotide database. Useful for gene discovery, sequence annotation, and comparative genomics.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py blastn -query my_query.fasta -db my_database -out results.out\033[0m\n'
    '     This command aligns sequences from \033[1;33mmy_query.fasta\033[0m with \033[1;33mm\033[0m and saves the results to \033[1;33mresults.out\033[0m.\n\n'
    '  \033[1;32m2. makeblastdb:\033[0m Create BLAST databases\n'
    '     Converts raw sequence files into a searchable database format for BLAST. Supports various input types, including nucleotides and proteins.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py makeblastdb -in my_sequences.fasta -dbtype nucl -title "My Nucleotide Database" -out my_database\033[0m\n'
    '     This command creates a nucleotide database from \033[1;33mm\033[0m with the title \033[1;33m"My Nucleotide Database"\033[0m, and outputs it to \033[1;33mm\033[0m.\n\n'
    '  \033[1;32m3. blastp:\033[0m Perform protein-protein alignments\n'
    '     Aligns protein sequences against a protein database to find similar proteins and infer their functional and evolutionary relationships.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py blastp -query my_protein_query.fasta -db my_protein_database -out results.out\033[0m\n'
    '     This command aligns sequences from \033[1;33mm\033[0m with \033[1;33mm\033[0m and saves the results to \033[1;33mresults.out\033[0m.\n\n'
    '  \033[1;32m4. blastx:\033[0m Translate nucleotide sequences and perform protein-protein alignments\n'
    '     Translates nucleotide sequences into proteins and then performs protein-protein alignments against a protein database.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py blastx -query my_nucleotide_query.fasta -db my_protein_database -out results.out\033[0m\n'
    '     This command translates sequences from \033[1;33mm\033[0m and aligns them with \033[1;33mm\033[0m, saving the results to \033[1;33mresults.out\033[0m.\n\n'
    '  \033[1;32m5. tblastn:\033[0m Perform protein-nucleotide alignments\n'
    '     Aligns protein sequences against a nucleotide database by translating the nucleotide sequences into proteins.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py tblastn -query my_protein_query.fasta -db my_nucleotide_database -out results.out\033[0m\n'
    '     This command aligns protein sequences from \033[1;33mm\033[0m with the translated nucleotide database \033[1;33mm\033[0m, saving results to \033[1;33mresults.out\033[0m.\n\n'
    '  \033[1;32m6. tblastx:\033[0m Perform nucleotide-nucleotide alignments with translation\n'
    '     Translates both nucleotide sequences in the query and database into proteins and performs protein-protein alignments.\n'
    '     \033[1;36mExample usage:\033[0m\n'
    '       \033[1;33mpython pyBlast.py tblastx -query my_nucleotide_query.fasta -db my_nucleotide_database -out results.out\033[0m\n'
    '     This command translates both sequences and aligns them, saving the results to \033[1;33mresults.out\033[0m.\n\n'
    'Basic usage examples:\n'
    '  - \033[1;32mNucleotide alignments:\033[0m\n'
    '    \033[1;33mpython pyBlast.py blastn -query my_query.fasta -db my_database\033[0m\n'
    '  - \033[1;32mCreate a BLAST database:\033[0m\n'
    '    \033[1;33mpython pyBlast.py makeblastdb -in my_sequences.fasta -dbtype nucl\033[0m\n'
    '  - \033[1;32mProtein alignments:\033[0m\n'
    '    \033[1;33mpython pyBlast.py blastp -query my_protein_query.fasta -db my_protein_database\033[0m\n'
    '  - \033[1;32mProtein alignments with translated nucleotide sequences:\033[0m\n'
    '    \033[1;33mpython pyBlast.py blastx -query my_nucleotide_query.fasta -db my_protein_database\033[0m\n'
    '  - \033[1;32mProtein-nucleotide alignments:\033[0m\n'
    '    \033[1;33mpython pyBlast.py tblastn -query my_protein_query.fasta -db my_nucleotide_database\033[0m\n'
    '  - \033[1;32mTranslated nucleotide alignments:\033[0m\n'
    '    \033[1;33mpython pyBlast.py tblastx -query my_nucleotide_query.fasta -db my_nucleotide_database\033[0m\n\n'
    'For detailed help on each subcommand, use:\n'
    '  - \033[1;33mpython pyBlast.py blastn -h\033[0m\n'
    '  - \033[1;33mpython pyBlast.py makeblastdb -h\033[0m\n'
    '  - \033[1;33mpython pyBlast.py blastp -h\033[0m\n'
    '  - \033[1;33mpython pyBlast.py blastx -h\033[0m\n'
    '  - \033[1;33mpython pyBlast.py tblastn -h\033[0m\n'
    '  - \033[1;33mpython pyBlast.py tblastx -h\033[0m\n'
)


    def define_blastn_arguments(self):
        blastn_parser = self.subparsers.add_parser(
            'blastn',
            help='Perform nucleotide-nucleotide sequence alignments',
            formatter_class=argparse.RawTextHelpFormatter,
            description=(
                'The \033[1;32mblastn\033[0m command aligns nucleotide sequences against a nucleotide database to find similarities. '
                'It is particularly useful for identifying homologous sequences, discovering genes, and annotating genome data.\n\n'
                'Key Features:\n'
                '  - High-speed searches for nucleotide sequences against nucleotide databases.\n'
                '  - Various filtering and scoring options to refine search results.\n'
                '  - Output formats for easy integration with other tools and analysis pipelines.\n\n'
                'Typical Usage:\n'
                '  \033[1;33mpython pyBlast.py blastn -query query.fasta -db database -out results.txt -evalue 1e-10\033[0m\n'
                '  This command performs a BLAST search using sequences from \033[1;33mquery.fasta\033[0m against the database \033[1;33mdatabase\033[0m, '
                'with results saved to \033[1;33mresults.txt\033[0m and an E-value threshold of \033[1;33m1e-10\033[0m.\n\n'
                'Additional Options:\n'
                '  \033[1;33m-query\033[0m: Input file containing nucleotide sequences to search.\n'
                '  \033[1;33m-db\033[0m: Name of the nucleotide database to search against.\n'
                '  \033[1;33m-out\033[0m: File to write the results to.\n'
                '  \033[1;33m-evalue\033[0m: E-value threshold for reporting matches.\n'
                '  \033[1;33m-task\033[0m: Specific BLASTN task to use (e.g., megablast, dc-megablast).\n'
            )
        )

        blastn_parser.add_argument(
            '-help', 
            action='store_true', 
            help='\033[1;36mPrint USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\033[0m'
        )
        blastn_parser.add_argument(
            '-version', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPrint version number; ignore other arguments\033[0m'
        )
        blastn_parser.add_argument(
            '-query', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mInput file name\033[0m\nSpecify the file containing the query sequence.'
        )
        blastn_parser.add_argument(
            '-query_loc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLocation on the query sequence in 1-based offsets (Format: start-stop)\033[0m\nDefine the range on the query sequence to use.'
        )
        blastn_parser.add_argument(
            '-strand', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mQuery strand(s) to search against database/subject\033[0m\nSpecify which strand(s) of the query sequence to use in the search.'
        )
        blastn_parser.add_argument(
            '-task', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTask to execute\033[0m\nChoose the specific BLASTN task to perform, such as “blastn” or “megablast”.'
        )
        blastn_parser.add_argument(
            '-db', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBLAST database name\033[0m\nSpecify the name of the BLAST database to search against.'
        )
        blastn_parser.add_argument(
            '-out', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mOutput file name\033[0m\nDefine the file where the search results will be saved.'
        )
        blastn_parser.add_argument(
            '-evalue', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mExpectation value (E) threshold for saving hits\033[0m\nSet the E-value threshold for reporting hits. Lower values yield more stringent criteria.'
        )
        blastn_parser.add_argument(
            '-word_size', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mWord size for wordfinder algorithm (length of best perfect match)\033[0m\nSpecify the length of the word size for initial alignment hits.'
        )
        blastn_parser.add_argument(
            '-gapopen', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to open a gap\033[0m\nDefine the penalty for opening a gap in the alignment.'
        )
        blastn_parser.add_argument(
            '-gapextend', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to extend a gap\033[0m\nSpecify the penalty for extending an existing gap in the alignment.'
        )
        blastn_parser.add_argument(
            '-penalty', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPenalty for a nucleotide mismatch\033[0m\nSet the penalty for mismatching nucleotides in the alignment.'
        )
        blastn_parser.add_argument(
            '-reward', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mReward for a nucleotide match\033[0m\nSpecify the reward for matching nucleotides in the alignment.'
        )
        blastn_parser.add_argument(
            '-use_index', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mUse MegaBLAST database index\033[0m\nIndicate whether to use the MegaBLAST index for faster searches.'
        )
        blastn_parser.add_argument(
            '-index_name', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMegaBLAST database index name (deprecated; use only for old style indices)\033[0m\nSpecify the index name for MegaBLAST databases. This is deprecated in favor of newer methods.'
        )
        blastn_parser.add_argument(
            '-subject', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSubject sequence(s) to search\033[0m\nProvide the subject sequences for the BLAST search.'
        )
        blastn_parser.add_argument(
            '-subject_loc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLocation on the subject sequence in 1-based offsets (Format: start-stop)\033[0m\nSpecify the location on the subject sequence for the search.'
        )
        blastn_parser.add_argument(
            '-outfmt', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36malignment view options\033[0m\nDefine the format for displaying alignments, such as tabular, XML, etc.'
        )
        blastn_parser.add_argument(
            '-show_gis', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mShow NCBI GIs in deflines?\033[0m\nSpecify whether to include NCBI GI numbers in the output definition lines.'
        )
        blastn_parser.add_argument(
            '-num_descriptions', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of database sequences to show one-line descriptions for\033[0m\nSet the number of sequences to display one-line descriptions for in the output.'
        )
        blastn_parser.add_argument(
            '-num_alignments', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of database sequences to show alignments for\033[0m\nSpecify the number of sequences to show alignments for in the output.'
        )
        blastn_parser.add_argument(
            '-line_length', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLine length for formatting alignments\033[0m\nDefine the length of lines for alignment formatting.'
        )
        blastn_parser.add_argument(
            '-html', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mProduce HTML output?\033[0m\nIndicate whether to generate HTML-formatted output for easier viewing in a web browser.'
        )
        blastn_parser.add_argument(
            '-sorthits', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSorting option for hits\033[0m\nSpecify how to sort the hits in the output, such as by score or E-value.'
        )
        blastn_parser.add_argument(
            '-sorthsps', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSorting option for HSPs\033[0m\nChoose the sorting option for High-Scoring Pairs (HSPs) in the output.'
        )
        blastn_parser.add_argument(
            '-dust', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFilter query sequence with DUST (Format: "yes", "level window linker", or "no" to disable)\033[0m\nSpecify whether to use DUST to filter the query sequence for low-complexity regions.'
        )
        blastn_parser.add_argument(
            '-filtering_db', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBLAST database containing filtering elements (i.e.: repeats)\033[0m\nProvide a database with filtering elements to use during the search.'
        )
        blastn_parser.add_argument(
            '-window_masker_taxid', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEnable WindowMasker filtering using a Taxonomic ID\033[0m\nSpecify a taxonomic ID for filtering with WindowMasker.'
        )
        blastn_parser.add_argument(
            '-window_masker_db', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEnable WindowMasker filtering using this repeats database.\033[0m\nProvide a database for WindowMasker filtering.'
        )
        blastn_parser.add_argument(
            '-soft_masking', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mApply filtering locations as soft masks\033[0m\nChoose whether to apply masking locations as soft masks.'
        )
        blastn_parser.add_argument(
            '-lcase_masking', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mUse lower case filtering in query and subject sequence(s)?\033[0m\nSpecify whether to apply lower case masking to query and subject sequences.'
        )
        blastn_parser.add_argument(
            '-gilist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to list of GIs\033[0m\nLimit the search to the specified list of GenInfo Identifiers (GIs).'
        )
        blastn_parser.add_argument(
            '-seqidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to list of SeqIDs\033[0m\nLimit the search to the specified list of sequence IDs (SeqIDs).'
        )
        blastn_parser.add_argument(
            '-negative_gilist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to everything except the specified GIs\033[0m\nExclude the specified list of GenInfo Identifiers (GIs) from the search.'
        )
        blastn_parser.add_argument(
            '-negative_seqidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to everything except the specified SeqIDs\033[0m\nExclude the specified list of sequence IDs (SeqIDs) from the search.'
        )
        blastn_parser.add_argument(
            '-taxids', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to include only the specified taxonomy IDs and their descendants (multiple IDs delimited by ",")\033[0m\nLimit the search to the specified taxonomy IDs and their descendants.'
        )
        blastn_parser.add_argument(
            '-negative_taxids', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to everything except the specified taxonomy IDs and their descendants (multiple IDs delimited by ",")\033[0m\nExclude the specified taxonomy IDs and their descendants from the search.'
        )
        blastn_parser.add_argument(
            '-taxidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to include only the specified taxonomy IDs and their descendants\033[0m\nLimit the search to the specified list of taxonomy IDs and their descendants.'
        )
        blastn_parser.add_argument(
            '-negative_taxidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search of database to everything except the specified taxonomy IDs and their descendants\033[0m\nExclude the specified list of taxonomy IDs and their descendants from the search.'
        )
        blastn_parser.add_argument(
            '-no_taxid_expansion', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mDo not expand the taxonomy IDs provided to their descendant taxonomy IDs\033[0m\nSpecify whether to disable expansion of taxonomy IDs to include their descendants.'
        )
        blastn_parser.add_argument(
            '-entrez_query', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search with the given Entrez query\033[0m\nApply an Entrez query to restrict the search results.'
        )
        blastn_parser.add_argument(
            '-db_soft_mask', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFiltering algorithm ID to apply to the BLAST database as soft masking\033[0m\nSpecify a filtering algorithm for soft masking of the BLAST database.'
        )
        blastn_parser.add_argument(
            '-db_hard_mask', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFiltering algorithm ID to apply to the BLAST database as hard masking\033[0m\nSpecify a filtering algorithm for hard masking of the BLAST database.'
        )
        blastn_parser.add_argument(
            '-perc_identity', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPercent identity\033[0m\nSet the minimum percent identity for reporting alignments.'
        )
        blastn_parser.add_argument(
            '-qcov_hsp_perc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPercent query coverage per HSP\033[0m\nSpecify the percent coverage of the query sequence required for High-Scoring Pairs (HSPs).'
        )
        blastn_parser.add_argument(
            '-max_hsps', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSet maximum number of HSPs per subject sequence to save for each query\033[0m\nLimit the number of HSPs per subject sequence to retain for each query.'
        )
        blastn_parser.add_argument(
            '-culling_limit', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mIf the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit\033[0m\nSpecify the culling limit for removing lower-scoring hits enveloped by higher-scoring ones.'
        )
        blastn_parser.add_argument(
            '-best_hit_overhang', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBest Hit algorithm overhang value (recommended value: 0.1)\033[0m\nSet the overhang value for the Best Hit algorithm, with a recommended value of 0.1.'
        )
        blastn_parser.add_argument(
            '-best_hit_score_edge', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBest Hit algorithm score edge value (recommended value: 0.1)\033[0m\nDefine the score edge value for the Best Hit algorithm, with a recommended value of 0.1.'
        )
        blastn_parser.add_argument(
            '-subject_besthit', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTurn on best hit per subject sequence\033[0m\nEnable the Best Hit algorithm for each subject sequence.'
        )
        blastn_parser.add_argument(
            '-max_target_seqs', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMaximum number of aligned sequences to keep\033[0m\nSet the maximum number of aligned sequences to retain in the output.'
        )
        blastn_parser.add_argument(
            '-template_type', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mDiscontiguous MegaBLAST template type\033[0m\nSpecify the template type for Discontiguous MegaBLAST.'
        )
        blastn_parser.add_argument(
            '-template_length', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mDiscontiguous MegaBLAST template length\033[0m\nSet the template length for Discontiguous MegaBLAST.'
        )
        blastn_parser.add_argument(
            '-dbsize', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEffective length of the database\033[0m\nDefine the effective length of the BLAST database.'
        )
        blastn_parser.add_argument(
            '-searchsp', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEffective length of the search space\033[0m\nSpecify the effective length of the search space for the BLAST search.'
        )
        blastn_parser.add_argument(
            '-import_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSearch strategy to use\033[0m\nProvide a search strategy for the BLAST search.'
        )
        blastn_parser.add_argument(
            '-export_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFile name to record the search strategy used\033[0m\nSpecify a file to record the search strategy employed.'
        )
        blastn_parser.add_argument(
            '-xdrop_ungap', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for ungapped extensions\033[0m\nSet the X-dropoff value in bits for ungapped extensions.'
        )
        blastn_parser.add_argument(
            '-xdrop_gap', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for preliminary gapped extensions\033[0m\nSpecify the X-dropoff value in bits for preliminary gapped extensions.'
        )
        blastn_parser.add_argument(
            '-xdrop_gap_final', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for final gapped alignment\033[0m\nDefine the X-dropoff value in bits for final gapped alignment.'
        )
        blastn_parser.add_argument(
            '-no_greedy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mUse non-greedy dynamic programming extension\033[0m\nSpecify whether to use a non-greedy approach for dynamic programming extension.'
        )
        blastn_parser.add_argument(
            '-min_raw_gapped_score', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMinimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages\033[0m\nSet the minimum raw gapped score required to retain an alignment.'
        )
        blastn_parser.add_argument(
            '-ungapped', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPerform ungapped alignment only?\033[0m\nSpecify whether to perform only ungapped alignment.'
        )
        blastn_parser.add_argument(
            '-window_size', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMultiple hits window size, use 0 to specify 1-hit algorithm\033[0m\nSet the window size for multiple hits; use 0 for a 1-hit algorithm.'
        )
        blastn_parser.add_argument(
            '-off_diagonal_range', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of off-diagonals to search for the 2nd hit, use 0 to turn off\033[0m\nSpecify the number of off-diagonals to search for the second hit; use 0 to turn off.'
        )
        blastn_parser.add_argument(
            '-parse_deflines', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mShould the query and subject defline(s) be parsed?\033[0m\nIndicate whether to parse the query and subject defline(s).'
        )
        blastn_parser.add_argument(
            '-num_threads', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of threads (CPUs) to use in the BLAST search\033[0m\nSpecify the number of threads or CPUs to use for the BLAST search.'
        )
        blastn_parser.add_argument(
            '-mt_mode', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMulti-thread mode to use in BLAST search\033[0m\nChoose the multi-thread mode for the BLAST search.'
        )
        blastn_parser.add_argument(
            '-remote', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mExecute search remotely?\033[0m\nSpecify whether to execute the BLAST search remotely.'
        )


    def define_makeblastdb_arguments(self):
        makeblastdb_parser = self.subparsers.add_parser(
            'makeblastdb',
            help='Create a BLAST database from input sequence files',
            formatter_class=argparse.RawTextHelpFormatter,
            description=(
                'The \033[1;32mmakeblastdb\033[0m command is used to generate a BLAST database from a given input sequence file. '
                'This database can then be used for sequence alignment searches using various BLAST tools. \n\n'
                'Key Features:\n'
                '  - Converts raw sequence data into a searchable database format.\n'
                '  - Supports nucleotide and protein sequences.\n'
                '  - Allows for custom database titles and descriptions.\n\n'
                'Typical Usage:\n'
                '  \033[1;33mpython pyBlast.py makeblastdb -in sequences.fasta -dbtype nucl -out my_database\033[0m\n'
                '  This command reads sequences from \033[1;33msequences.fasta\033[0m and creates a nucleotide database named \033[1;33mmy_database\033[0m.\n\n'
                'Additional Options:\n'
                '  \033[1;33m-in\033[0m: Input file containing the sequences.\n'
                '  \033[1;33m-dbtype\033[0m: Type of database (nucl for nucleotide, prot for protein).\n'
                '  \033[1;33m-out\033[0m: Output database file name.\n'
                '  \033[1;33m-title\033[0m: Custom title for the database.\n'
                '  \033[1;33m-taxid\033[0m: Taxonomy ID for the sequences, if applicable.\n'
            )
        )
        makeblastdb_parser.add_argument(
            '-help', 
            action='store_true', 
            help='\033[1;33mPrint USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\033[0m'
        )
        makeblastdb_parser.add_argument(
            '-in', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mInput file/database name\033[0m\nSpecify the file or database name that you want to use as input.'
        )
        makeblastdb_parser.add_argument(
            '-input_type', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mType of the data specified in input_file\033[0m\nDefine the type of data in the input file, such as nucleotide or protein sequences.'
        )
        makeblastdb_parser.add_argument(
            '-dbtype', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMolecule type of target db\033[0m\nSpecify the type of molecule for the target database (e.g., nucl or prot).'
        )
        makeblastdb_parser.add_argument(
            '-title', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTitle for BLAST database\033[0m\nProvide a title or description for the BLAST database being created.'
        )
        makeblastdb_parser.add_argument(
            '-parse_seqids', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mOption to parse seqid for FASTA input\033[0m\nIf enabled for FASTA input, seqids are parsed manually; otherwise, they are parsed automatically for other input types.'
        )
        makeblastdb_parser.add_argument(
            '-hash_index', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCreate index of sequence hash values\033[0m\nGenerate a hash index to speed up sequence lookup operations.'
        )
        makeblastdb_parser.add_argument(
            '-mask_data', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mComma-separated list of input files containing masking data\033[0m\nProvide files with masking data, produced by NCBI masking applications, to be used during database creation.'
        )
        makeblastdb_parser.add_argument(
            '-mask_id', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mComma-separated list of strings to uniquely identify the masking algorithm\033[0m\nList identifiers for different masking algorithms applied to the data.'
        )
        makeblastdb_parser.add_argument(
            '-mask_desc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mComma-separated list of free form strings to describe the masking algorithm details\033[0m\nProvide descriptions for the masking algorithms used, which can be free text.'
        )
        makeblastdb_parser.add_argument(
            '-gi_mask', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCreate GI indexed masking data\033[0m\nGenerate masking data indexed by GI numbers (GenInfo Identifier).'
        )
        makeblastdb_parser.add_argument(
            '-gi_mask_name', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mComma-separated list of masking data output files\033[0m\nSpecify the output files for GI indexed masking data.'
        )
        makeblastdb_parser.add_argument(
            '-out', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mName of BLAST database to be created\033[0m\nName the output BLAST database file to be created.'
        )
        makeblastdb_parser.add_argument(
            '-blastdb_version', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mVersion of BLAST database to be created\033[0m\nSpecify the version of the BLAST database format to use (e.g., 5 or 6).'
        )
        makeblastdb_parser.add_argument(
            '-max_file_sz', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMaximum file size for BLAST database files\033[0m\nDefine the maximum file size for each BLAST database file.'
        )
        makeblastdb_parser.add_argument(
            '-metadata_output_prefix', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPath prefix for location of database files in metadata\033[0m\nProvide the directory path where metadata and other database files will be stored.'
        )
        makeblastdb_parser.add_argument(
            '-logfile', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFile to which the program log should be redirected\033[0m\nSpecify a file to capture the log output for troubleshooting and analysis.'
        )
        makeblastdb_parser.add_argument(
            '-taxid', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTaxonomy ID to assign to all sequences\033[0m\nAssign a specific taxonomy ID to all sequences in the database.'
        )
        makeblastdb_parser.add_argument(
            '-taxid_map', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mText file mapping sequence IDs to taxonomy IDs\033[0m\nProvide a file that maps sequence IDs to their corresponding taxonomy IDs for correct annotation.'
        )
        makeblastdb_parser.add_argument(
            '-oid_masks', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36m0x01 Exclude Model\033[0m\nApply specific model exclusion settings to the database creation process.'
        )


    def define_blastp_arguments(self):
        blastp_parser = self.subparsers.add_parser(
            'blastp',
            help='Perform protein-protein sequence alignments',
            formatter_class=argparse.RawTextHelpFormatter,
            description=(
                'The \033[1;32mblastp\033[0m command performs protein-protein alignments to identify similar protein sequences. '
                'It compares an input protein sequence against a protein database to find homologous proteins and infer their functional and evolutionary relationships.\n\n'
                'Key Features:\n'
                '  - Searches for protein matches in a protein database.\n'
                '  - Provides detailed alignment results including percent identity and alignment scores.\n'
                '  - Flexible output formats for integration into various analysis workflows.\n\n'
                'Typical Usage:\n'
                '  \033[1;33mpython pyBlast.py blastp -query protein_query.fasta -db protein_db -out results.txt -evalue 1e-5\033[0m\n'
                '  This command aligns sequences from \033[1;33mprotein_query.fasta\033[0m against \033[1;33mprotein_db\033[0m, '
                'saving results to \033[1;33mresults.txt\033[0m with an E-value threshold of \033[1;33m1e-5\033[0m.\n\n'
                'Additional Options:\n'
                '  \033[1;33m-query\033[0m: File containing protein sequences to align.\n'
                '  \033[1;33m-db\033[0m: Name of the protein database to search against.\n'
                '  \033[1;33m-out\033[0m: File to write the results to.\n'
                '  \033[1;33m-evalue\033[0m: E-value threshold for reporting matches.\n'
                '  \033[1;33m-matrix\033[0m: Scoring matrix to use for alignment (e.g., BLOSUM62).\n'
            )
        )
        blastp_parser.add_argument(
            '-help', 
            action='store_true', 
            help='\033[1;33mPrint USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\033[0m'
        )
        blastp_parser.add_argument(
            '-version', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;33mPrint version number; ignore other arguments\033[0m\nDisplays the current version of the `blastp` tool.'
        )
        blastp_parser.add_argument(
            '-query', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mInput file name\033[0m\nSpecify the input file containing protein sequences to query against the database.'
        )
        blastp_parser.add_argument(
            '-query_loc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLocation on the query sequence in 1-based offsets (Format: start-stop)\033[0m\nDefine the region of the query sequence to be used, in the format start-stop.'
        )
        blastp_parser.add_argument(
            '-subject', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSubject sequence(s) to search\033[0m\nSpecify the subject sequence or sequences for the search.'
        )
        blastp_parser.add_argument(
            '-subject_loc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLocation on the subject sequence in 1-based offsets (Format: start-stop)\033[0m\nDefine the region of the subject sequence to be used, in the format start-stop.'
        )
        blastp_parser.add_argument(
            '-task', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTask to execute\033[0m\nSpecify the BLASTP task to be executed, such as “blastp” or other custom tasks.'
        )
        blastp_parser.add_argument(
            '-db', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBLAST database name\033[0m\nSpecify the name of the BLAST database to search against.'
        )
        blastp_parser.add_argument(
            '-dbsize', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEffective length of the database\033[0m\nSet the effective length of the database, which affects the E-value calculations.'
        )
        blastp_parser.add_argument(
            '-evalue', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mExpectation value (E) threshold for saving hits\033[0m\nSpecify the E-value threshold for reporting hits, where lower values denote more significant matches.'
        )
        blastp_parser.add_argument(
            '-word_size', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mWord size for wordfinder algorithm (length of best perfect match)\033[0m\nSet the length of the words used to find initial matches during the search.'
        )
        blastp_parser.add_argument(
            '-gapopen', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to open a gap\033[0m\nSpecify the penalty for opening a gap in the alignment.'
        )
        blastp_parser.add_argument(
            '-gapextend', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to extend a gap\033[0m\nSpecify the penalty for extending a gap in the alignment.'
        )
        blastp_parser.add_argument(
            '-qcov_hsp_perc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPercent query coverage per HSP\033[0m\nDefine the minimum percent coverage of the query sequence for each High-Scoring Pair (HSP) to be reported.'
        )
        blastp_parser.add_argument(
            '-max_hsps', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSet maximum number of HSPs per subject sequence to save for each query\033[0m\nLimit the number of HSPs per subject sequence to save for each query.'
        )
        blastp_parser.add_argument(
            '-xdrop_ungap', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for ungapped extensions\033[0m\nSet the X-dropoff value in bits for ungapped extensions during the alignment process.'
        )
        blastp_parser.add_argument(
            '-xdrop_gap', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for preliminary gapped extensions\033[0m\nSpecify the X-dropoff value in bits for preliminary gapped extensions.'
        )
        blastp_parser.add_argument(
            '-xdrop_gap_final', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mX-dropoff value (in bits) for final gapped alignment\033[0m\nSet the X-dropoff value in bits for the final gapped alignment.'
        )
        blastp_parser.add_argument(
            '-searchsp', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEffective length of the search space\033[0m\nSpecify the effective length of the search space for the BLAST search.'
        )
        blastp_parser.add_argument(
            '-seg', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSEG options for masking\033[0m\nProvide options for masking regions of low complexity in the sequences using the SEG algorithm.'
        )
        blastp_parser.add_argument(
            '-soft_masking', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mApply filtering locations as soft masks\033[0m\nSpecify whether to apply soft masking to the filtered locations in the query and subject sequences.'
        )
        blastp_parser.add_argument(
            '-matrix', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mScoring matrix for protein sequences\033[0m\nSpecify the scoring matrix to use for protein sequence alignments (e.g., BLOSUM62, PAM250).'
        )
        blastp_parser.add_argument(
            '-threshold', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mThreshold for reporting a match\033[0m\nSet the threshold score required for reporting a match in the BLAST results.'
        )
        blastp_parser.add_argument(
            '-culling_limit', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCulling limit for hit filtering\033[0m\nDefine the culling limit for filtering out less significant hits based on the alignment score.'
        )
        blastp_parser.add_argument(
            '-best_hit_overhang', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBest Hit algorithm overhang value (recommended value: 0.1)\033[0m\nSpecify the overhang value for the Best Hit algorithm to control the threshold for reporting the best hits.'
        )
        blastp_parser.add_argument(
            '-best_hit_score_edge', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBest Hit algorithm score edge value (recommended value: 0.1)\033[0m\nDefine the score edge value for the Best Hit algorithm to determine the cutoff score for reporting hits.'
        )
        blastp_parser.add_argument(
            '-subject_besthit', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTurn on best hit per subject sequence\033[0m\nEnable the option to report only the best hit per subject sequence.'
        )
        blastp_parser.add_argument(
            '-window_size', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMultiple hits window size, use 0 to specify 1-hit algorithm\033[0m\nSet the window size for multiple hits. Use 0 for a 1-hit algorithm.'
        )
        blastp_parser.add_argument(
            '-lcase_masking', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mUse lower case filtering in query and subject sequence(s)?\033[0m\nSpecify whether to apply lower case masking to the query and subject sequences.'
        )
        blastp_parser.add_argument(
            '-parse_deflines', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mShould the query and subject defline(s) be parsed?\033[0m\nDetermine whether to parse the definition lines of the query and subject sequences.'
        )
        blastp_parser.add_argument(
            '-outfmt', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mAlignment view options\033[0m\nSpecify the format for the alignment output, such as tabular, XML, or others.'
        )
        blastp_parser.add_argument(
            '-show_gis', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mShow NCBI GIs in deflines?\033[0m\nIndicate whether to include NCBI GI numbers in the output definition lines.'
        )
        blastp_parser.add_argument(
            '-num_descriptions', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of database sequences to show one-line descriptions for\033[0m\nSet the number of database sequences to display one-line descriptions for in the output.'
        )
        blastp_parser.add_argument(
            '-num_alignments', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of database sequences to show alignments for\033[0m\nSpecify the number of database sequences to display alignments for in the output.'
        )
        blastp_parser.add_argument(
            '-line_length', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLine length for formatting alignments\033[0m\nDefine the length of lines for formatting alignment output.'
        )
        blastp_parser.add_argument(
            '-html', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mProduce HTML output?\033[0m\nSpecify whether to generate HTML-formatted output for easier viewing in a web browser.'
        )
        blastp_parser.add_argument(
            '-sorthits', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSorting option for hits\033[0m\nChoose the sorting order for the hits in the output, such as by score or E-value.'
        )
        blastp_parser.add_argument(
            '-sorthsps', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSorting option for HSPs\033[0m\nSpecify the sorting option for High-Scoring Pairs (HSPs) in the output.'
        )
        blastp_parser.add_argument(
            '-max_target_seqs', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMaximum number of aligned sequences to keep\033[0m\nSet the maximum number of aligned sequences to retain in the output.'
        )
        blastp_parser.add_argument(
            '-num_threads', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mNumber of threads (CPUs) to use in the BLAST search\033[0m\nSpecify the number of threads or CPUs to use for parallel processing during the BLAST search.'
        )
        blastp_parser.add_argument(
            '-mt_mode', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mMulti-thread mode to use in BLAST search\033[0m\nChoose the multi-threading mode for the BLAST search to optimize performance.'
        )
        blastp_parser.add_argument(
            '-ungapped', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPerform ungapped alignment only?\033[0m\nSpecify whether to perform only ungapped alignments without introducing gaps.'
        )
        blastp_parser.add_argument(
            '-remote', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mExecute search remotely?\033[0m\nIndicate whether to perform the search remotely on a BLAST server, rather than locally.'
        )
        blastp_parser.add_argument(
            '-import_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSearch strategy to use\033[0m\nSpecify a predefined search strategy to be used during the search.'
        )
        blastp_parser.add_argument(
            '-export_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFile name to record the search strategy used\033[0m\nProvide a file name to save the search strategy used in the current run.'
        )
        blastp_parser.add_argument(
            '-comp_based_stats', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mComposition-based statistics\033[0m\nEnable composition-based statistics to improve sensitivity in detecting weak matches.'
        )
        blastp_parser.add_argument(
            '-use_sw_tback', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mUse Smith-Waterman traceback\033[0m\nSpecify whether to use the Smith-Waterman algorithm for traceback to improve alignment accuracy.'
        )


    def define_blastx_arguments(self):
        blastx_parser = self.subparsers.add_parser(
            'blastx',
            help='Perform protein query vs. translated nucleotide database search',
            formatter_class=argparse.RawTextHelpFormatter,
            description=(
                'The \033[1;32mblastx\033[0m command aligns protein queries against a nucleotide database that is dynamically translated '
                'in all reading frames, allowing for comparisons across species and genomic analysis. This can be particularly useful for discovering protein coding genes and annotating transcriptomic data.\n\n'
                'Key Features:\n'
                '  - Translate nucleotide sequences dynamically for protein-level comparisons.\n'
                '  - Comprehensive filtering and scoring mechanisms to refine your results.\n'
                '  - Customizable output formats for various downstream applications and analyses.\n\n'
                'Typical Usage:\n'
                '  \033[1;33mpython pyBlast.py blastx -query query.fasta -db protein_db -out results.txt -evalue 1e-5\033[0m\n'
                '  This command performs a BLASTX search using sequences from \033[1;33mquery.fasta\033[0m against the protein database \033[1;33mprotein_db\033[0m, '
                'with results saved to \033[1;33mresults.txt\033[0m and an E-value threshold of \033[1;33m1e-5\033[0m.\n\n'
                'Additional Options:\n'
                '  \033[1;33m-query\033[0m: Input file containing nucleotide sequences to be translated and searched.\n'
                '  \033[1;33m-db\033[0m: Name of the protein database to search against.\n'
                '  \033[1;33m-out\033[0m: File to write the results to.\n'
                '  \033[1;33m-evalue\033[0m: E-value threshold for reporting matches.\n'
                '  \033[1;33m-task\033[0m: Specific BLASTX task to use (e.g., blastx-fast, blastx).\n'
            )
        )

        blastx_parser.add_argument(
            '-help', 
            action='store_true', 
            help='\033[1;36mPrint USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\033[0m\nDisplay the usage message, program description, and argument information, ignoring any other specified parameters.'
        )
        blastx_parser.add_argument(
            '-version', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mPrint version number; ignore other arguments\033[0m\nShow the version number of the program and ignore all other parameters.'
        )
        blastx_parser.add_argument(
            '-import_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSearch strategy file to import\033[0m\nSpecify a search strategy file to import for the BLAST search.'
        )
        blastx_parser.add_argument(
            '-export_search_strategy', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFile to record the search strategy used\033[0m\nProvide a file name to export and save the search strategy used in the BLAST search.'
        )
        blastx_parser.add_argument(
            '-task', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mTask to execute\033[0m\nSpecify the BLAST task to perform, such as blastx, blastx-fast, etc.'
        )
        blastx_parser.add_argument(
            '-db', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mBLAST database name\033[0m\nProvide the name of the BLAST database to search against.'
        )
        blastx_parser.add_argument(
            '-dbsize', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mEffective length of the database\033[0m\nSpecify the effective length of the database, used for statistical calculations.'
        )
        blastx_parser.add_argument(
            '-gilist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to list of GIs\033[0m\nLimit the search to a specified list of GI (GenInfo Identifier) numbers.'
        )
        blastx_parser.add_argument(
            '-seqidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to list of SeqIDs\033[0m\nLimit the search to a specified list of sequence IDs.'
        )
        blastx_parser.add_argument(
            '-negative_gilist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to everything except the specified GIs\033[0m\nExclude a list of GI numbers from the search.'
        )
        blastx_parser.add_argument(
            '-negative_seqidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to everything except the specified SeqIDs\033[0m\nExclude a list of sequence IDs from the search.'
        )
        blastx_parser.add_argument(
            '-taxids', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to specified taxonomy IDs and descendants\033[0m\nLimit the search to specific taxonomy IDs and their descendant taxa.'
        )
        blastx_parser.add_argument(
            '-negative_taxids', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to everything except the specified taxonomy IDs and descendants\033[0m\nExclude specified taxonomy IDs and their descendant taxa from the search.'
        )
        blastx_parser.add_argument(
            '-taxidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to list of taxonomy IDs and descendants\033[0m\nLimit the search to a list of taxonomy IDs and their descendants.'
        )
        blastx_parser.add_argument(
            '-negative_taxidlist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to everything except the specified taxonomy IDs and descendants\033[0m\nExclude a list of taxonomy IDs and their descendants from the search.'
        )
        blastx_parser.add_argument(
            '-no_taxid_expansion', 
            action='store_true', 
            help='\033[1;36mDo not expand the taxonomy IDs provided to their descendant taxonomy IDs\033[0m\nPrevent the automatic expansion of taxonomy IDs to include their descendants.'
        )
        blastx_parser.add_argument(
            '-ipglist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to list of IPGs\033[0m\nLimit the search to a specified list of IPGs (Identical Protein Groups).'
        )
        blastx_parser.add_argument(
            '-negative_ipglist', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search to everything except the specified IPGs\033[0m\nExclude a list of IPGs from the search.'
        )
        blastx_parser.add_argument(
            '-entrez_query', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mRestrict search with the given Entrez query\033[0m\nFilter the search results using an Entrez query.'
        )
        blastx_parser.add_argument(
            '-db_soft_mask', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFiltering algorithm to apply to the database as soft masking\033[0m\nApply soft masking to the database using the specified filtering algorithm.'
        )
        blastx_parser.add_argument(
            '-db_hard_mask', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mFiltering algorithm to apply to the database as hard masking\033[0m\nApply hard masking to the database using the specified filtering algorithm.'
        )
        blastx_parser.add_argument(
            '-subject', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mSubject sequence(s) to search\033[0m\nProvide the subject sequence(s) for the BLAST search.'
        )
        blastx_parser.add_argument(
            '-subject_loc', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mLocation on the subject sequence in 1-based offsets (Format: start-stop)\033[0m\nSpecify the location on the subject sequence using 1-based offsets (e.g., 100-200).'
        )
        blastx_parser.add_argument(
            '-query', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mInput file name\033[0m\nSpecify the input file name containing the query sequence(s).'
        )
        blastx_parser.add_argument(
            '-out', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mOutput file name\033[0m\nProvide the output file name for saving the results of the BLAST search.'
        )
        blastx_parser.add_argument(
            '-evalue', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mExpectation value (E) threshold for saving hits\033[0m\nSet the expectation value (E-value) threshold for retaining significant hits.'
        )
        blastx_parser.add_argument(
            '-word_size', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mWord size for wordfinder algorithm\033[0m\nSpecify the word size for the wordfinder algorithm in the BLAST search.'
        )
        blastx_parser.add_argument(
            '-gapopen', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to open a gap\033[0m\nDefine the penalty for opening a gap in alignments.'
        )
        blastx_parser.add_argument(
            '-gapextend', 
            nargs='?', 
            default='NotUsed', 
            const='UseWithoutParam', 
            help='\033[1;36mCost to extend a gap\033[0m\nDefine the penalty for extending a gap in alignments.'
        )


        def define_tblastn_arguments(self):
            tblastn_parser = self.subparsers.add_parser(
                'tblastn',
                help='Perform protein query against translated nucleotide database search',
                formatter_class=argparse.RawTextHelpFormatter,
                description=(
                    'The \033[1;32mtblastn\033[0m command aligns protein sequences against a nucleotide database that is dynamically translated into all reading frames. '
                    'This is particularly useful for detecting potential coding regions in nucleotide sequences and for comparative genomics.\n\n'
                    'Key Features:\n'
                    '  - Translates nucleotide databases into protein sequences dynamically for alignment against protein queries.\n'
                    '  - Various filtering, scoring, and alignment customization options to refine search results.\n'
                    '  - Supports different output formats for downstream analysis and integration with other bioinformatics tools.\n\n'
                    'Typical Usage:\n'
                    '  \033[1;33mpython pyBlast.py tblastn -query protein.fasta -db nucleotide_db -out results.txt -evalue 1e-3\033[0m\n'
                    '  This command performs a TBLASTN search using sequences from \033[1;33mprotein.fasta\033[0m against the nucleotide database \033[1;33mnucleotide_db\033[0m, '
                    'with results saved to \033[1;33mresults.txt\033[0m and an E-value threshold of \033[1;33m1e-3\033[0m.\n\n'
                    'Additional Options:\n'
                    '  \033[1;33m-query\033[0m: Input file containing protein sequences to search against the translated nucleotide database.\n'
                    '  \033[1;33m-db\033[0m: Name of the nucleotide database to search, which will be dynamically translated into protein sequences.\n'
                    '  \033[1;33m-out\033[0m: File to write the results to.\n'
                    '  \033[1;33m-evalue\033[0m: E-value threshold for reporting significant matches.\n'
                    '  \033[1;33m-task\033[0m: Specific TBLASTN task to execute (e.g., tblastn, tblastn-fast).\n'
                )
            )

            # Adding tblastn specific arguments
            tblastn_parser.add_argument(
                '-help',
                action='store_true',
                help=(
                    '\033[1;33m-help\033[0m\n'
                    'Display help information for the tblastn command, including usage, description, and options.'
                )
            )

            tblastn_parser.add_argument(
                '-version',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-version\033[0m\n'
                    'Print the version number of the TBLASTN tool. This option will ignore all other parameters and only print the version information.'
                )
            )

            tblastn_parser.add_argument(
                '-import_search_strategy',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-import_search_strategy\033[0m\n'
                    'Specify a file from which to import search strategy parameters.'
                )
            )

            tblastn_parser.add_argument(
                '-export_search_strategy',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-export_search_strategy\033[0m\n'
                    'Specify a file to export the search strategy parameters used in the search.'
                )
            )

            tblastn_parser.add_argument(
                '-task',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-task\033[0m\n'
                    'Specify the TBLASTN task to execute, such as `tblastn` for standard searches or `tblastn-fast` for faster searches.'
                )
            )

            tblastn_parser.add_argument(
                '-db',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db\033[0m\n'
                    'Specify the nucleotide database to search against. The database will be dynamically translated into protein sequences for the search.'
                )
            )

            tblastn_parser.add_argument(
                '-dbsize',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-dbsize\033[0m\n'
                    'Specify the effective length of the database for calculations.'
                )
            )

            tblastn_parser.add_argument(
                '-gilist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-gilist\033[0m\n'
                    'Restrict the search to a specific list of GI numbers provided in a file.'
                )
            )

            tblastn_parser.add_argument(
                '-seqidlist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-seqidlist\033[0m\n'
                    'Restrict the search to a specific list of SeqIDs provided in a file.'
                )
            )

            tblastn_parser.add_argument(
                '-negative_gilist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-negative_gilist\033[0m\n'
                    'Restrict the search to everything except the GI numbers specified in the list file.'
                )
            )

            tblastn_parser.add_argument(
                '-negative_seqidlist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-negative_seqidlist\033[0m\n'
                    'Restrict the search to everything except the SeqIDs specified in the list file.'
                )
            )

            tblastn_parser.add_argument(
                '-taxids',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-taxids\033[0m\n'
                    'Restrict the search to specified taxonomy IDs and their descendants.'
                )
            )

            tblastn_parser.add_argument(
                '-negative_taxids',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-negative_taxids\033[0m\n'
                    'Restrict the search to everything except the specified taxonomy IDs and their descendants.'
                )
            )

            tblastn_parser.add_argument(
                '-taxidlist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-taxidlist\033[0m\n'
                    'Restrict the search to a list of taxonomy IDs and their descendants provided in a file.'
                )
            )

            tblastn_parser.add_argument(
                '-negative_taxidlist',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-negative_taxidlist\033[0m\n'
                    'Restrict the search to everything except the taxonomy IDs and their descendants listed in the file.'
                )
            )

            tblastn_parser.add_argument(
                '-no_taxid_expansion',
                action='store_true',
                help=(
                    '\033[1;33m-no_taxid_expansion\033[0m\n'
                    'Do not expand the provided taxonomy IDs to include their descendant taxonomy IDs.'
                )
            )

            tblastn_parser.add_argument(
                '-entrez_query',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-entrez_query\033[0m\n'
                    'Restrict the search with the given Entrez query.'
                )
            )

            tblastn_parser.add_argument(
                '-db_soft_mask',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db_soft_mask\033[0m\n'
                    'Specify a filtering algorithm to apply to the database as soft masking.'
                )
            )

            tblastn_parser.add_argument(
                '-db_hard_mask',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db_hard_mask\033[0m\n'
                    'Specify a filtering algorithm to apply to the database as hard masking.'
                )
            )

            tblastn_parser.add_argument(
                '-subject',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-subject\033[0m\n'
                    'Specify the input file containing subject sequences for the search.'
                )
            )

            tblastn_parser.add_argument(
                '-subject_loc',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-subject_loc\033[0m\n'
                    'Specify the location on the subject sequence in 1-based offsets (Format: start-stop).'
                )
            )

            tblastn_parser.add_argument(
                '-query',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-query\033[0m\n'
                    'Specify the input file containing protein sequences to search against the translated nucleotide database.'
                )
            )

            tblastn_parser.add_argument(
                '-out',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-out\033[0m\n'
                    'Specify the file to write the results of the search to.'
                )
            )

            tblastn_parser.add_argument(
                '-evalue',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-evalue\033[0m\n'
                    'Specify the Expectation value (E) threshold for saving hits.'
                )
            )

            tblastn_parser.add_argument(
                '-word_size',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-word_size\033[0m\n'
                    'Specify the word size for the wordfinder algorithm.'
                )
            )

            tblastn_parser.add_argument(
                '-gapopen',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-gapopen\033[0m\n'
                    'Specify the cost to open a gap.'
                )
            )

            tblastn_parser.add_argument(
                '-gapextend',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-gapextend\033[0m\n'
                    'Specify the cost to extend a gap.'
                )
            )

            tblastn_parser.add_argument(
                '-qcov_hsp_perc',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-qcov_hsp_perc\033[0m\n'
                    'Specify the percent query coverage per HSP.'
                )
            )

            tblastn_parser.add_argument(
                '-max_hsps',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_hsps\033[0m\n'
                    'Specify the maximum number of HSPs per subject sequence to save for each query.'
                )
            )

            tblastn_parser.add_argument(
                '-xdrop_ungap',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-xdrop_ungap\033[0m\n'
                    'Specify the X-dropoff value for ungapped extensions.'
                )
            )

            tblastn_parser.add_argument(
                '-xdrop_gap',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-xdrop_gap\033[0m\n'
                    'Specify the X-dropoff value for preliminary gapped extensions.'
                )
            )

            tblastn_parser.add_argument(
                '-xdrop_gap_final',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-xdrop_gap_final\033[0m\n'
                    'Specify the X-dropoff value for final gapped alignment.'
                )
            )

            tblastn_parser.add_argument(
                '-searchsp',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-searchsp\033[0m\n'
                    'Specify the effective length of the search space.'
                )
            )

            tblastn_parser.add_argument(
                '-sum_stats',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sum_stats\033[0m\n'
                    'Enable or disable sum statistics.'
                )
            )

            tblastn_parser.add_argument(
                '-db_gencode',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db_gencode\033[0m\n'
                    'Specify the genetic code to use for database translation.'
                )
            )

            tblastn_parser.add_argument(
                '-ungapped',
                action='store_true',
                help=(
                    '\033[1;33m-ungapped\033[0m\n'
                    'Perform ungapped alignment only.'
                )
            )

            tblastn_parser.add_argument(
                '-max_intron_length',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_intron_length\033[0m\n'
                    'Specify the maximum intron length.'
                )
            )

            tblastn_parser.add_argument(
                '-seg',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-seg\033[0m\n'
                    'Specify SEG filtering options.'
                )
            )

            tblastn_parser.add_argument(
                '-soft_masking',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-soft_masking\033[0m\n'
                    'Apply filtering locations as soft masks.'
                )
            )

            tblastn_parser.add_argument(
                '-matrix',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-matrix\033[0m\n'
                    'Specify the scoring matrix name.'
                )
            )

            tblastn_parser.add_argument(
                '-threshold',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-threshold\033[0m\n'
                    'Specify the minimum score threshold for extending hits.'
                )
            )

            tblastn_parser.add_argument(
                '-culling_limit',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-culling_limit\033[0m\n'
                    'Specify the limit on hits enveloped by higher-scoring hits.'
                )
            )

            tblastn_parser.add_argument(
                '-best_hit_overhang',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-best_hit_overhang\033[0m\n'
                    'Specify the best hit algorithm overhang value.'
                )
            )

            tblastn_parser.add_argument(
                '-best_hit_score_edge',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-best_hit_score_edge\033[0m\n'
                    'Specify the best hit algorithm score edge value.'
                )
            )

            tblastn_parser.add_argument(
                '-subject_besthit',
                action='store_true',
                help=(
                    '\033[1;33m-subject_besthit\033[0m\n'
                    'Enable the best hit per subject sequence.'
                )
            )

            tblastn_parser.add_argument(
                '-window_size',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-window_size\033[0m\n'
                    'Specify the window size for multiple hits (use 0 to specify 1-hit algorithm).'
                )
            )

            tblastn_parser.add_argument(
                '-lcase_masking',
                action='store_true',
                help=(
                    '\033[1;33m-lcase_masking\033[0m\n'
                    'Enable lowercase masking in query and subject sequences.'
                )
            )

            tblastn_parser.add_argument(
                '-query_loc',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-query_loc\033[0m\n'
                    'Specify the location on the query sequence in 1-based offsets (Format: start-stop).'
                )
            )

            tblastn_parser.add_argument(
                '-parse_deflines',
                action='store_true',
                help=(
                    '\033[1;33m-parse_deflines\033[0m\n'
                    'Enable parsing of query and subject deflines.'
                )
            )

            tblastn_parser.add_argument(
                '-outfmt',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-outfmt\033[0m\n'
                    'Specify the output alignment format.'
                )
            )

            tblastn_parser.add_argument(
                '-show_gis',
                action='store_true',
                help=(
                    '\033[1;33m-show_gis\033[0m\n'
                    'Display NCBI GI numbers in deflines.'
                )
            )

            tblastn_parser.add_argument(
                '-num_descriptions',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_descriptions\033[0m\n'
                    'Specify the number of descriptions to show in the output.'
                )
            )

            tblastn_parser.add_argument(
                '-num_alignments',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_alignments\033[0m\n'
                    'Specify the number of alignments to show in the output.'
                )
            )

            tblastn_parser.add_argument(
                '-line_length',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-line_length\033[0m\n'
                    'Specify the line length for the output.'
                )
            )

            tblastn_parser.add_argument(
                '-html',
                action='store_true',
                help=(
                    '\033[1;33m-html\033[0m\n'
                    'Generate output in HTML format.'
                )
            )

            tblastn_parser.add_argument(
                '-sorthits',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sorthits\033[0m\n'
                    'Specify how to sort hits in the output.'
                )
            )

            tblastn_parser.add_argument(
                '-sorthsps',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sorthsps\033[0m\n'
                    'Specify how to sort HSPs in the output.'
                )
            )

            tblastn_parser.add_argument(
                '-max_target_seqs',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_target_seqs\033[0m\n'
                    'Specify the maximum number of target sequences to report.'
                )
            )

            tblastn_parser.add_argument(
                '-num_threads',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_threads\033[0m\n'
                    'Specify the number of threads to use for parallel processing.'
                )
            )

            tblastn_parser.add_argument(
                '-mt_mode',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-mt_mode\033[0m\n'
                    'Specify the mode of multi-threading for BLAST searches.'
                )
            )

            tblastn_parser.add_argument(
                '-remote',
                action='store_true',
                help=(
                    '\033[1;33m-remote\033[0m\n'
                    'Perform the search using the NCBI BLAST service instead of local resources.'
                )
            )

            tblastn_parser.add_argument(
                '-comp_based_stats',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-comp_based_stats\033[0m\n'
                    'Enable or disable composition-based statistics.'
                )
            )

            tblastn_parser.add_argument(
                '-use_sw_tback',
                action='store_true',
                help=(
                    '\033[1;33m-use_sw_tback\033[0m\n'
                    'Use Smith-Waterman traceback for alignments.'
                )
            )

            tblastn_parser.add_argument(
                '-in_pssm',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-in_pssm\033[0m\n'
                    'Specify a PSI-BLAST checkpoint file for use in the search.'
                )
            )


        def define_tblastx_arguments(self):
            tblastx_parser = self.subparsers.add_parser(
                'tblastx',
                help='Perform translated nucleotide-nucleotide sequence alignments',
                formatter_class=argparse.RawTextHelpFormatter,
                description=(
                    'The \033[1;32mtblastx\033[0m command aligns nucleotide sequences against a nucleotide database by translating both query and database sequences into proteins. '
                    'It is useful for detecting homologous genes and annotating genomes when both sequences are nucleotide-based.\n\n'
                    'Key Features:\n'
                    '  - Translated nucleotide-to-nucleotide alignments.\n'
                    '  - Various filtering and scoring options to enhance search results.\n'
                    '  - Output formats for integration with other tools and pipelines.\n\n'
                    'Typical Usage:\n'
                    '  \033[1;33mpython pyBlast.py tblastx -query query.fasta -db database -out results.txt -evalue 1e-10\033[0m\n'
                    '  This command performs a BLAST search using translated sequences from \033[1;33mquery.fasta\033[0m against the translated nucleotide database \033[1;33mdatabase\033[0m, '
                    'with results saved to \033[1;33mresults.txt\033[0m and an E-value threshold of \033[1;33m1e-10\033[0m.\n\n'
                    'Additional Options:\n'
                    '  \033[1;33m-db\033[0m: Name of the nucleotide database to search against.\n'
                    '  \033[1;33m-query\033[0m: Input file containing nucleotide sequences to translate and search.\n'
                    '  \033[1;33m-out\033[0m: File to write the search results.\n'
                    '  \033[1;33m-evalue\033[0m: E-value threshold for reporting matches.\n'
                    '  \033[1;33m-word_size\033[0m: Word size for the wordfinder algorithm.\n'
                    '  \033[1;33m-qcov_hsp_perc\033[0m: Percent query coverage per HSP.\n'
                    '  \033[1;33m-max_hsps\033[0m: Maximum number of HSPs per subject sequence to save for each query.\n'
                    '  \033[1;33m-xdrop_ungap\033[0m: X-dropoff value for ungapped extensions.\n'
                    '  \033[1;33m-searchsp\033[0m: Effective length of the search space.\n'
                    '  \033[1;33m-sum_stats\033[0m: Enable or disable sum statistics.\n'
                    '  \033[1;33m-max_intron_length\033[0m: Maximum intron length.\n'
                    '  \033[1;33m-seg\033[0m: SEG filtering options.\n'
                    '  \033[1;33m-soft_masking\033[0m: Apply filtering locations as soft masks.\n'
                    '  \033[1;33m-matrix\033[0m: Scoring matrix name.\n'
                    '  \033[1;33m-threshold\033[0m: Minimum score threshold for extending hits.\n'
                    '  \033[1;33m-culling_limit\033[0m: Limit on hits enveloped by higher-scoring hits.\n'
                    '  \033[1;33m-best_hit_overhang\033[0m: Best hit algorithm overhang value.\n'
                    '  \033[1;33m-best_hit_score_edge\033[0m: Best hit algorithm score edge value.\n'
                    '  \033[1;33m-subject_besthit\033[0m: Enable the best hit per subject sequence.\n'
                    '  \033[1;33m-window_size\033[0m: Window size for multiple hits (use 0 for 1-hit algorithm).\n'
                    '  \033[1;33m-lcase_masking\033[0m: Enable lowercase masking in query and subject sequences.\n'
                    '  \033[1;33m-query_loc\033[0m: Location on the query sequence in 1-based offsets (Format: start-stop).\n'
                    '  \033[1;33m-strand\033[0m: Strand to search against (e.g., plus, minus).\n'
                    '  \033[1;33m-parse_deflines\033[0m: Enable parsing of query and subject deflines.\n'
                    '  \033[1;33m-query_gencode\033[0m: Genetic code to use for query sequence translation.\n'
                    '  \033[1;33m-db_gencode\033[0m: Genetic code to use for database sequence translation.\n'
                    '  \033[1;33m-outfmt\033[0m: Output alignment format.\n'
                    '  \033[1;33m-show_gis\033[0m: Display NCBI GI numbers in deflines.\n'
                    '  \033[1;33m-num_descriptions\033[0m: Number of descriptions to show in the output.\n'
                    '  \033[1;33m-num_alignments\033[0m: Number of alignments to show in the output.\n'
                    '  \033[1;33m-line_length\033[0m: Line length for the output.\n'
                    '  \033[1;33m-html\033[0m: Generate output in HTML format.\n'
                    '  \033[1;33m-sorthits\033[0m: How to sort hits in the output.\n'
                    '  \033[1;33m-sorthsps\033[0m: How to sort HSPs in the output.\n'
                    '  \033[1;33m-max_target_seqs\033[0m: Maximum number of target sequences to report.\n'
                    '  \033[1;33m-num_threads\033[0m: Number of threads to use for parallel processing.\n'
                    '  \033[1;33m-remote\033[0m: Perform the search using the NCBI BLAST service instead of local resources.\n'
                )
            )

            tblastx_parser.add_argument(
                '-db',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db\033[0m\n'
                    'Specify the nucleotide database to search against.'
                )
            )

            tblastx_parser.add_argument(
                '-query',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-query\033[0m\n'
                    'Specify the input file containing nucleotide sequences to translate and search.'
                )
            )

            tblastx_parser.add_argument(
                '-out',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-out\033[0m\n'
                    'Specify the file to write the results of the search to.'
                )
            )

            tblastx_parser.add_argument(
                '-evalue',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-evalue\033[0m\n'
                    'Specify the Expectation value (E) threshold for saving hits.'
                )
            )

            tblastx_parser.add_argument(
                '-word_size',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-word_size\033[0m\n'
                    'Specify the word size for the wordfinder algorithm.'
                )
            )

            tblastx_parser.add_argument(
                '-qcov_hsp_perc',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-qcov_hsp_perc\033[0m\n'
                    'Specify the percent query coverage per HSP.'
                )
            )

            tblastx_parser.add_argument(
                '-max_hsps',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_hsps\033[0m\n'
                    'Specify the maximum number of HSPs per subject sequence to save for each query.'
                )
            )

            tblastx_parser.add_argument(
                '-xdrop_ungap',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-xdrop_ungap\033[0m\n'
                    'Specify the X-dropoff value for ungapped extensions.'
                )
            )

            tblastx_parser.add_argument(
                '-searchsp',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-searchsp\033[0m\n'
                    'Specify the effective length of the search space.'
                )
            )

            tblastx_parser.add_argument(
                '-sum_stats',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sum_stats\033[0m\n'
                    'Enable or disable sum statistics in the results.'
                )
            )

            tblastx_parser.add_argument(
                '-max_intron_length',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_intron_length\033[0m\n'
                    'Specify the maximum intron length to consider.'
                )
            )

            tblastx_parser.add_argument(
                '-seg',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-seg\033[0m\n'
                    'Specify SEG filtering options.'
                )
            )

            tblastx_parser.add_argument(
                '-soft_masking',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-soft_masking\033[0m\n'
                    'Apply filtering locations as soft masks.'
                )
            )

            tblastx_parser.add_argument(
                '-matrix',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-matrix\033[0m\n'
                    'Specify the scoring matrix to use.'
                )
            )

            tblastx_parser.add_argument(
                '-threshold',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-threshold\033[0m\n'
                    'Specify the minimum score threshold for extending hits.'
                )
            )

            tblastx_parser.add_argument(
                '-culling_limit',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-culling_limit\033[0m\n'
                    'Specify the limit on hits enveloped by higher-scoring hits.'
                )
            )

            tblastx_parser.add_argument(
                '-best_hit_overhang',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-best_hit_overhang\033[0m\n'
                    'Specify the best hit algorithm overhang value.'
                )
            )

            tblastx_parser.add_argument(
                '-best_hit_score_edge',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-best_hit_score_edge\033[0m\n'
                    'Specify the best hit algorithm score edge value.'
                )
            )

            tblastx_parser.add_argument(
                '-subject_besthit',
                action='store_true',
                help=(
                    '\033[1;33m-subject_besthit\033[0m\n'
                    'Enable the best hit per subject sequence.'
                )
            )

            tblastx_parser.add_argument(
                '-window_size',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-window_size\033[0m\n'
                    'Specify the window size for multiple hits (use 0 for 1-hit algorithm).'
                )
            )

            tblastx_parser.add_argument(
                '-lcase_masking',
                action='store_true',
                help=(
                    '\033[1;33m-lcase_masking\033[0m\n'
                    'Enable lowercase masking in query and subject sequences.'
                )
            )

            tblastx_parser.add_argument(
                '-query_loc',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-query_loc\033[0m\n'
                    'Specify location on the query sequence in 1-based offsets (Format: start-stop).'
                )
            )

            tblastx_parser.add_argument(
                '-strand',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-strand\033[0m\n'
                    'Specify the strand to search against (e.g., plus, minus).'
                )
            )

            tblastx_parser.add_argument(
                '-parse_deflines',
                action='store_true',
                help=(
                    '\033[1;33m-parse_deflines\033[0m\n'
                    'Enable parsing of query and subject deflines.'
                )
            )

            tblastx_parser.add_argument(
                '-query_gencode',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-query_gencode\033[0m\n'
                    'Specify the genetic code to use for query sequence translation.'
                )
            )

            tblastx_parser.add_argument(
                '-db_gencode',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-db_gencode\033[0m\n'
                    'Specify the genetic code to use for database sequence translation.'
                )
            )

            tblastx_parser.add_argument(
                '-outfmt',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-outfmt\033[0m\n'
                    'Specify the output alignment format.'
                )
            )

            tblastx_parser.add_argument(
                '-show_gis',
                action='store_true',
                help=(
                    '\033[1;33m-show_gis\033[0m\n'
                    'Display NCBI GI numbers in deflines.'
                )
            )

            tblastx_parser.add_argument(
                '-num_descriptions',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_descriptions\033[0m\n'
                    'Specify the number of descriptions to show in the output.'
                )
            )

            tblastx_parser.add_argument(
                '-num_alignments',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_alignments\033[0m\n'
                    'Specify the number of alignments to show in the output.'
                )
            )

            tblastx_parser.add_argument(
                '-line_length',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-line_length\033[0m\n'
                    'Specify the line length for the output.'
                )
            )

            tblastx_parser.add_argument(
                '-html',
                action='store_true',
                help=(
                    '\033[1;33m-html\033[0m\n'
                    'Generate output in HTML format.'
                )
            )

            tblastx_parser.add_argument(
                '-sorthits',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sorthits\033[0m\n'
                    'Specify how to sort hits in the output.'
                )
            )

            tblastx_parser.add_argument(
                '-sorthsps',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-sorthsps\033[0m\n'
                    'Specify how to sort HSPs in the output.'
                )
            )

            tblastx_parser.add_argument(
                '-max_target_seqs',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-max_target_seqs\033[0m\n'
                    'Specify the maximum number of target sequences to report.'
                )
            )

            tblastx_parser.add_argument(
                '-num_threads',
                nargs='?',
                default='NotUsed',
                const='UseWithoutParam',
                help=(
                    '\033[1;33m-num_threads\033[0m\n'
                    'Specify the number of threads to use for parallel processing.'
                )
            )

            tblastx_parser.add_argument(
                '-remote',
                action='store_true',
                help=(
                    '\033[1;33m-remote\033[0m\n'
                    'Perform the search using the NCBI BLAST service instead of local resources.'
                )
            )

            tblastx_parser.add_argument(
                '-version',
                action='store_true',
                help=(
                    '\033[1;33m-version\033[0m\n'
                    'Show the version of the tblastx tool.'
                )
            )


    def run(self):
        args = self.parser.parse_args()
        if args.command is None:
            self.parser.print_help()
            sys.exit(1)

        program_path = None
        if args.command == 'blastn':
            program_path = shutil.which('blastn')
        elif args.command == 'makeblastdb':
            program_path = shutil.which('makeblastdb')
        elif args.command == 'blastp':
            program_path = shutil.which('blastp')
        elif args.command == 'blastx':
            program_path = shutil.which('blastx')
        elif args.command == 'tblastn':
            program_path = shutil.which('tblastn')
        elif args.command == 'tblastx':
            program_path = shutil.which('tblastx')
        else:
            print(f"Unknown command: {args.command}")
            sys.exit(1)

        if program_path is None:
            print(f"Error: {args.command} not found. Please make sure it is installed and in your PATH.")
            sys.exit(1)

        command = [program_path]

        if args.command == 'blastn':
            self.run_blastn(args, command)
        elif args.command == 'makeblastdb':
            self.run_makeblastdb(args, command)
        elif args.command == 'blastp':
            self.run_blastp(args, command)
        elif args.command == 'blastx':
            self.run_blastx(args, command)
        elif args.command == 'tblastn':
            self.run_tblastn(args, command)
        elif args.command == 'tblastx':
            self.run_tblastx(args, command)

    def run_blastn(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'blastn':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

    def run_makeblastdb(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'makeblastdb':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

    def run_blastp(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'blastp':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

    def run_blastx(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'blastx':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

    def run_tblastn(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'tblastn':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

    def run_tblastx(self, args, command):
        #print("#########var(args):##########",vars(args).items())
        for arg_name, value in vars(args).items():
            if arg_name == 'help' and value is True:
                self.parser.print_help()
                sys.exit()
            elif arg_name == 'help' and value is False:
                continue
            elif value == 'NotUsed':
                continue
            elif value == 'UseWithoutParam':
                command.append(f"-{arg_name}")
            elif arg_name == 'command' and value == 'tblastx':
                continue
            else:
                command.append(f"-{arg_name}")
                command.append(f"{value}")
        subprocess.run(command)
        #print(command)

if __name__ == "__main__":
    tool = BlastTool()
    tool.define_blastn_arguments()
    tool.define_makeblastdb_arguments()
    tool.define_blastp_arguments()
    tool.define_blastx_arguments()   
    tool.define_tblastn_arguments()  
    tool.define_tblastx_arguments()  
    tool.run()

