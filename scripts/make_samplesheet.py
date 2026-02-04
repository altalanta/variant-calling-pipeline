#!/usr/bin/env python3

"""
Generate samplesheet for variant calling pipeline
Scans directories for FASTQ files and creates properly formatted samplesheet
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def find_fastq_files(directory, pattern=None):
    """Find FASTQ files in directory"""
    fastq_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    fastq_files = []

    directory = Path(directory)

    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")

    for ext in fastq_extensions:
        if pattern:
            files = directory.rglob(f"*{pattern}*{ext}")
        else:
            files = directory.rglob(f"*{ext}")
        fastq_files.extend(files)

    return sorted(fastq_files)


def parse_sample_info(filename):
    """Extract sample information from filename"""
    # Common patterns for sample names
    patterns = [
        r"^([^_]+)_.*_R?([12]).*\.f(ast)?q(\.gz)?$",  # sample_lib_R1.fastq.gz
        r"^([^_]+).*[._]([12])\.f(ast)?q(\.gz)?$",  # sample.1.fastq.gz
        r"^([^_.]+).*_R?([12])_.*\.f(ast)?q(\.gz)?$",  # sample_R1_001.fastq.gz
    ]

    filename = Path(filename).name

    for pattern in patterns:
        match = re.match(pattern, filename)
        if match:
            sample_id = match.group(1)
            read_num = int(match.group(2))
            return sample_id, read_num

    # Fallback: try to extract sample name before first separator
    sample_id = filename.split("_")[0].split(".")[0]

    # Try to determine read number from filename
    if "_R1" in filename or ".1." in filename or "_1." in filename:
        read_num = 1
    elif "_R2" in filename or ".2." in filename or "_2." in filename:
        read_num = 2
    else:
        read_num = None

    return sample_id, read_num


def pair_fastq_files(fastq_files):
    """Group FASTQ files into pairs"""
    samples = defaultdict(dict)
    unpaired = []

    for fastq_file in fastq_files:
        sample_id, read_num = parse_sample_info(fastq_file)

        if read_num is None:
            unpaired.append(fastq_file)
            continue

        if read_num in samples[sample_id]:
            print(f"Warning: Multiple R{read_num} files found for sample {sample_id}")
            print(f"  Existing: {samples[sample_id][read_num]}")
            print(f"  New: {fastq_file}")
            continue

        samples[sample_id][read_num] = fastq_file

    # Check for complete pairs
    paired_samples = {}
    for sample_id, reads in samples.items():
        if 1 in reads and 2 in reads:
            paired_samples[sample_id] = {"read1": reads[1], "read2": reads[2]}
        else:
            missing = [str(i) for i in [1, 2] if i not in reads]
            print(f"Warning: Sample {sample_id} missing R{','.join(missing)} file(s)")
            if reads:
                unpaired.extend(reads.values())

    if unpaired:
        print(f"Warning: {len(unpaired)} unpaired files found:")
        for f in unpaired[:5]:  # Show first 5
            print(f"  {f}")
        if len(unpaired) > 5:
            print(f"  ... and {len(unpaired) - 5} more")

    return paired_samples


def detect_platform_and_library(fastq_file):
    """Detect sequencing platform and library from filename or path"""
    filename = str(fastq_file)

    # Platform detection
    platform = "ILLUMINA"  # Default
    if any(x in filename.upper() for x in ["NOVASEQ", "NEXTSEQ", "HISEQ", "MISEQ"]):
        platform = "ILLUMINA"
    elif "PACBIO" in filename.upper() or "PB" in filename.upper():
        platform = "PACBIO"
    elif "NANOPORE" in filename.upper() or "ONT" in filename.upper():
        platform = "ONT"

    # Library detection (try to extract from path or filename)
    parts = Path(fastq_file).parts
    library = None

    # Look for library indicators in path
    for part in parts:
        if "lib" in part.lower():
            library = part
            break
        elif re.match(r"L\d+", part):  # Lane pattern
            library = part
            break

    # Extract lane from filename
    lane_match = re.search(r"L(\d+)", filename)
    lane = lane_match.group(1) if lane_match else "1"

    return platform, library, lane


def create_samplesheet(paired_samples, output_file, auto_detect=True):
    """Create samplesheet CSV file"""

    with open(output_file, "w", newline="") as csvfile:
        fieldnames = ["sample_id", "read1", "read2", "platform", "library", "lane"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for sample_id, reads in sorted(paired_samples.items()):
            row = {
                "sample_id": sample_id,
                "read1": str(reads["read1"]),
                "read2": str(reads["read2"]),
            }

            if auto_detect:
                platform, library, lane = detect_platform_and_library(reads["read1"])
                row["platform"] = platform
                row["library"] = library or sample_id
                row["lane"] = lane
            else:
                row["platform"] = "ILLUMINA"
                row["library"] = sample_id
                row["lane"] = "1"

            writer.writerow(row)

    print(f"Samplesheet written to: {output_file}")
    print(f"Found {len(paired_samples)} paired samples")


def validate_samplesheet(samplesheet_file):
    """Validate the generated samplesheet"""
    required_columns = ["sample_id", "read1", "read2"]
    issues = []

    try:
        with open(samplesheet_file) as f:
            reader = csv.DictReader(f)

            # Check header
            if not all(col in reader.fieldnames for col in required_columns):
                missing = [col for col in required_columns if col not in reader.fieldnames]
                issues.append(f"Missing required columns: {missing}")
                return issues

            # Check rows
            for i, row in enumerate(reader, 1):
                # Check required fields
                for col in required_columns:
                    if not row[col].strip():
                        issues.append(f"Row {i}: Empty {col}")

                # Check file existence
                for fastq_col in ["read1", "read2"]:
                    fastq_path = row[fastq_col].strip()
                    if fastq_path and not Path(fastq_path).exists():
                        issues.append(f"Row {i}: File not found: {fastq_path}")

                # Check for duplicate sample IDs
                # (Would need to collect all rows first for this check)

    except Exception as e:
        issues.append(f"Error reading samplesheet: {e}")

    return issues


def main():
    parser = argparse.ArgumentParser(
        description="Generate samplesheet for variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan current directory for FASTQ files
  python make_samplesheet.py --input .

  # Scan specific directory with pattern
  python make_samplesheet.py --input /data/fastq --pattern "*WGS*"
  
  # Custom output file
  python make_samplesheet.py --input /data --output my_samples.csv
  
  # Disable auto-detection of platform/library
  python make_samplesheet.py --input /data --no-auto-detect
""",
    )

    parser.add_argument(
        "--input", "-i", required=True, help="Input directory containing FASTQ files"
    )

    parser.add_argument(
        "--output",
        "-o",
        default="samplesheet.csv",
        help="Output samplesheet file (default: samplesheet.csv)",
    )

    parser.add_argument("--pattern", "-p", help="Pattern to match in filenames (optional)")

    parser.add_argument(
        "--no-auto-detect",
        action="store_true",
        help="Disable auto-detection of platform and library info",
    )

    parser.add_argument(
        "--validate", action="store_true", help="Validate the generated samplesheet"
    )

    parser.add_argument(
        "--dry-run", action="store_true", help="Show what would be done without creating files"
    )

    args = parser.parse_args()

    try:
        # Find FASTQ files
        print(f"Scanning directory: {args.input}")
        if args.pattern:
            print(f"Using pattern: {args.pattern}")

        fastq_files = find_fastq_files(args.input, args.pattern)

        if not fastq_files:
            print("No FASTQ files found!")
            sys.exit(1)

        print(f"Found {len(fastq_files)} FASTQ files")

        # Pair files
        paired_samples = pair_fastq_files(fastq_files)

        if not paired_samples:
            print("No paired samples found!")
            sys.exit(1)

        if args.dry_run:
            print("\nDry run - would create samplesheet with:")
            for sample_id, reads in sorted(paired_samples.items()):
                print(f"  {sample_id}: {reads['read1']} + {reads['read2']}")
            sys.exit(0)

        # Create samplesheet
        create_samplesheet(paired_samples, args.output, auto_detect=not args.no_auto_detect)

        # Validate if requested
        if args.validate:
            print("\nValidating samplesheet...")
            issues = validate_samplesheet(args.output)
            if issues:
                print("Validation issues found:")
                for issue in issues:
                    print(f"  - {issue}")
                sys.exit(1)
            else:
                print("Samplesheet validation passed!")

        print(f"\nSamplesheet created: {args.output}")
        print("Review the file and edit if needed before running the pipeline.")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
