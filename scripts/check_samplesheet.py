#!/usr/bin/env python3

"""
Samplesheet validation script for variant calling pipeline
Validates input samplesheet format and checks file accessibility
"""

import argparse
import csv
import sys
from pathlib import Path


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample_id,read1,read2,platform,library,lane
    SAMPLE_1,SAMPLE_1_R1.fastq.gz,SAMPLE_1_R2.fastq.gz,ILLUMINA,LIB_1,1
    SAMPLE_2,SAMPLE_2_R1.fastq.gz,SAMPLE_2_R2.fastq.gz,ILLUMINA,LIB_2,2
    """

    required_columns = ["sample_id", "read1", "read2"]
    optional_columns = ["platform", "library", "lane"]
    valid_platforms = ["ILLUMINA", "PACBIO", "ONT", "IONTORRENT"]

    sample_mapping_dict = {}

    with open(file_in) as fin:
        ## Check header
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            print_error("Header not found or file is empty.", "Line", 1)

        # Normalize header (strip whitespace and convert to lowercase for comparison)
        normalized_header = [col.strip().lower() for col in reader.fieldnames]
        required_normalized = [col.lower() for col in required_columns]

        # Check for required columns
        missing_columns = [col for col in required_columns if col.lower() not in normalized_header]
        if missing_columns:
            print_error(f"Missing required column(s): {', '.join(missing_columns)}", "Header", "")

        # Check for unknown columns
        all_valid_columns = required_columns + optional_columns
        valid_normalized = [col.lower() for col in all_valid_columns]
        unknown_columns = [
            reader.fieldnames[i]
            for i, col in enumerate(normalized_header)
            if col not in valid_normalized
        ]
        if unknown_columns:
            print(
                f"Warning: Unknown column(s) found: {', '.join(unknown_columns)}", file=sys.stderr
            )

        ## Check sample entries
        sample_ids = set()
        for line_num, row in enumerate(reader, 2):  # Start at line 2 (after header)
            # Check for empty row
            if all(value.strip() == "" for value in row.values()):
                continue

            # Get sample_id (required)
            sample_id = row.get("sample_id", "").strip()
            if not sample_id:
                print_error("sample_id cannot be empty.", "Line", line_num)

            # Check for duplicate sample_ids
            if sample_id in sample_ids:
                print_error(f"Duplicate sample_id found: {sample_id}", "Line", line_num)
            sample_ids.add(sample_id)

            # Check read files (required)
            read1 = row.get("read1", "").strip()
            read2 = row.get("read2", "").strip()

            if not read1:
                print_error("read1 file path cannot be empty.", "Line", line_num)
            if not read2:
                print_error("read2 file path cannot be empty.", "Line", line_num)

            # Validate file extensions
            valid_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
            if read1 and not any(read1.endswith(ext) for ext in valid_extensions):
                print_error(
                    f"read1 file has invalid extension: {read1}. Expected: {', '.join(valid_extensions)}",
                    "Line",
                    line_num,
                )
            if read2 and not any(read2.endswith(ext) for ext in valid_extensions):
                print_error(
                    f"read2 file has invalid extension: {read2}. Expected: {', '.join(valid_extensions)}",
                    "Line",
                    line_num,
                )

            # Check file accessibility
            if read1:
                read1_path = Path(read1)
                if not read1_path.exists():
                    print_error(f"read1 file does not exist: {read1}", "Line", line_num)
                elif not read1_path.is_file():
                    print_error(f"read1 path is not a file: {read1}", "Line", line_num)

            if read2:
                read2_path = Path(read2)
                if not read2_path.exists():
                    print_error(f"read2 file does not exist: {read2}", "Line", line_num)
                elif not read2_path.is_file():
                    print_error(f"read2 path is not a file: {read2}", "Line", line_num)

            # Validate optional fields
            platform = row.get("platform", "ILLUMINA").strip()
            if platform and platform.upper() not in valid_platforms:
                print(
                    f"Warning: Unknown platform '{platform}' on line {line_num}. Valid platforms: {', '.join(valid_platforms)}",
                    file=sys.stderr,
                )
                platform = "ILLUMINA"  # Default to ILLUMINA

            library = row.get("library", sample_id).strip()
            if not library:
                library = sample_id  # Default to sample_id

            lane = row.get("lane", "1").strip()
            if not lane:
                lane = "1"  # Default to lane 1
            elif not lane.isdigit():
                print(
                    f"Warning: Lane '{lane}' is not numeric on line {line_num}. Setting to '1'.",
                    file=sys.stderr,
                )
                lane = "1"

            ## Create sample mapping dictionary
            sample_mapping_dict[sample_id] = {
                "read1": read1,
                "read2": read2,
                "platform": platform.upper(),
                "library": library,
                "lane": lane,
            }

    ## Write validated samplesheet
    if sample_mapping_dict:
        with open(file_out, "w") as fout:
            writer = csv.writer(fout)
            # Write header
            writer.writerow(["sample_id", "read1", "read2", "platform", "library", "lane"])

            # Write sample data
            for sample_id, sample_data in sorted(sample_mapping_dict.items()):
                writer.writerow(
                    [
                        sample_id,
                        sample_data["read1"],
                        sample_data["read2"],
                        sample_data["platform"],
                        sample_data["library"],
                        sample_data["lane"],
                    ]
                )

        print(f"Validated samplesheet written to: {file_out}", file=sys.stderr)
        print(f"Found {len(sample_mapping_dict)} valid samples", file=sys.stderr)
    else:
        print_error("No valid samples found in samplesheet", "File", file_in)


def print_error(message, label, value):
    """Print error message and exit"""
    print(f"ERROR: {message}", file=sys.stderr)
    if label and value:
        print(f"{label}: {value}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Validate samplesheet for variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Required samplesheet format:
    sample_id,read1,read2
    NA12878,/path/to/NA12878_R1.fastq.gz,/path/to/NA12878_R2.fastq.gz

Optional columns (with defaults):
    platform (default: ILLUMINA)
    library (default: sample_id)  
    lane (default: 1)

Examples:
    # Basic validation
    python check_samplesheet.py input.csv output.csv
    
    # Full format with optional columns
    sample_id,read1,read2,platform,library,lane
    NA12878,NA12878_R1.fq.gz,NA12878_R2.fq.gz,ILLUMINA,WGS_lib1,1
""",
    )

    parser.add_argument("file_in", metavar="FILE_IN", help="Input samplesheet CSV file")

    parser.add_argument(
        "file_out", metavar="FILE_OUT", help="Output validated samplesheet CSV file"
    )

    parser.add_argument(
        "--check-files",
        action="store_true",
        help="Check if FASTQ files exist and are accessible (default: enabled)",
    )

    parser.add_argument(
        "--no-check-files",
        action="store_true",
        help="Skip checking if FASTQ files exist (for testing)",
    )

    args = parser.parse_args()

    # Validate input file exists
    if not Path(args.file_in).exists():
        print_error(f"Input samplesheet not found: {args.file_in}", "", "")

    # Check if input file is readable
    try:
        with open(args.file_in) as f:
            pass
    except Exception as e:
        print_error(f"Cannot read input samplesheet: {e}", "", "")

    # Create output directory if needed
    output_path = Path(args.file_out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        check_samplesheet(args.file_in, args.file_out)
        print("Samplesheet validation completed successfully!", file=sys.stderr)
    except Exception as e:
        print_error(f"Samplesheet validation failed: {e}", "", "")


if __name__ == "__main__":
    main()
