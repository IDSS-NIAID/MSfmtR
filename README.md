# MSfmtR


See test files for a example workflows (need to include an example
here…). The configuration details are below.

## Configuration

### PD Filtering Defaults

- Directory parameters
  - `input_dir`: character, path to directory containing raw data
    (default: ‘.’).
  - `fasta_dir`: character path to directory containing fasta files
    (default: ‘.’).
  - `output_dir`: character, path to directory for output (default:
    ‘.’).
- Input file parameters
  - `in_file`: character, name of raw data file (default:
    `list.files(input_dir, pattern = 'xlsx')`),
  - `in_sheet`: character, name of the sheet in the raw data file (imports the first sheet by default),
- Sample parameters
  - `conditions`: character, vector of the conditions to be used in the
    analysis (default: NULL). If `conditions` is null, an attempt will be made
    to guess the conditions from the data.

### Spectronaut DIA Defaults

Any configuration options not defined in the yaml file will be set to
the package defaults. The default values are:

- Directory parameters
  - `input_dir`: character, path to directory containing raw data
    (default: ‘.’).
  - `fasta_dir`: character path to directory containing fasta files
    (default: ‘.’),
  - `output_dir`: character, path to directory for output (default:
    ‘.’),
- Input file parameters
  - `in_file`: character, name of raw data file (default:
    `list.files(input_dir, pattern = 'tsv')`),
  - `in_delim`: character, delimiter for raw data file (default: '\t')
- Output file parameters
  - `out_xlsx`: character, name of excel output file (default:
    `gsub('.tsv', '.xlsx', in_file, fixed = TRUE)`).
  - `out_sqlite`: character, name of SQLite output file (default:
    `gsub('.tsv', '.sqlite', in_file, fixed = TRUE)`).
  - `sheet`: character, name of the excel sheet for output to `out_xlsx` (default:
    `gsub('.tsv', '', in_file, fixed = TRUE)`).
- Sample parameters
  - `smp_grp_rep`: character, column name in raw data file specifying
    each sample / group / replicate (default: ‘R.FileName’).
  - `merge_method`: character, method used to summarize / group
    abundances among samples / groups / replicates (default: ‘median’).
- Filtering parameters
  - `uloq`: numeric, upper limit of quantification (default: Inf).
  - `lloq`: numeric, lower limit of quantification (default: 0).
  - `cont_fasta`: character path to fasta file containing contaminants
    (default:
    `system.file('extdata/Universal_Contaminants.fasta', package = 'MSfmtR')`),
  - `max_ratio`: numeric, maximum threshold to report. Any ratio above
    `max_ratio` or below `1 / max_ratio` will be truncated at this level
    (default: 100).
- Metadata parameters
  - `fasta_meta`: character, path to fasta metadata (default:
    `list.files(fasta_dir, pattern = 'fasta') |> grep(pattern = cont_fasta, invert = TRUE, value = TRUE)`).
  - `taxId`: character, default taxonomy ID to use for protein accession
    ID lookup. (default: 9606).
- MSStats parameters
  - `ratios`: character, vector of the contrasts to be used in the
    MSstats analysis (default: NULL). If ratios is null, all conditions
    will be compared based on the `R.Condition` column in `in_file`.
    Each ratio in the list should be of the form “<group1>/<group2>”, so
    for group1=Case and group2=Control, the ratio would be
    “Case/Control”.
  - `groups`: character, vector of the groups from `R.Condition` to be
    used in the creation of `ratios` (default: NULL). This can be used
    to specify a subgroup of conditions to compare (all conditions are
    used by default). If `groups` is defined and `ratios` is null, all
    combinations of `groups` will be used to construct a list of ratios
    to calculate.
  - `peptide_summary`: character flagging level of peptide summarization
    (ignored when `preprocess` is TRUE). Options: ‘PEP’ for peptide /
    elution group level (default), ‘FG’ for fragment group level, or
    ‘none’ to leave at the fragment level.
  - `normMeasure`: character, normalization measure (ignored when
    `preprocess` is TRUE or when `peptide_summary != 'none'`, default:
    ‘NormalizedPeakArea’).
  - `preprocess`: logical, preprocess data using MSstats (default:
    FALSE).
  - `format`: character, format of data (see `?raw_to_fld` for options,
    default: ‘MSstats’).
- Formatting parameters
  - `protein_header_fill`: character, color for protein header rows in
    excel output (default: “\#A7CDF0”).
  - `protein_rows_fill`: character, color for prtein rows in excel
    output (default: “\#DDEBF7”).
  - `peptide_header_fill`: character, color for peptide header rows in
    excel output (default: “\#F0CBA8”)
  - `peptide_rows_fill`: character, color for peptide rows in excel
    output (default: “\#FCE4D6”).
- Checkpoint parameters
  - `checkpoints`: character vector, specifying which intermediate files
    should be saved for debug purposes (default:
    `c('xlsx', 'sql', 'processed', 'protein', 'peptide', 'wb')`).
  - `processed_checkpoint`: character name the intermediate checkpoint
    file for processed data (default:
    `paste0(config$sheet, '_processed.RData')`),
  - `peptide_checkpoint`: character, name of the intermeditate
    checkpoint file for peptide data (default:
    `paste0(config$sheet, '_peptide.RData')`).
  - `protein_checkpoint`: character, name of the intermediate checkpoint
    file for protein data (default:
    `paste0(config$sheet, '_protein.RData')`).
  - `wb_checkpoint`: character, name of the intermedaiate checkpoint
    file for workbook data (default:
    `paste0(config$sheet, '_wb.RData')`).
