#!/usr/bin/env python3

"""
Verify pipeline installation and dependencies
Quick check to ensure all components are properly set up
"""

import argparse
import subprocess
import sys
from pathlib import Path


def check_command(command, version_flag="--version"):
    """Check if a command is available and get version"""
    try:
        result = subprocess.run([command, version_flag], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version = result.stdout.strip().split("\n")[0]
            return True, version
        else:
            return False, result.stderr.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError) as e:
        return False, str(e)


def check_file_exists(filepath, description):
    """Check if a file exists"""
    path = Path(filepath)
    exists = path.exists()
    size = path.stat().st_size if exists else 0
    return exists, f"Size: {size} bytes" if exists else "File not found"


def check_docker():
    """Check Docker availability"""
    # Check if Docker is installed
    docker_available, docker_version = check_command("docker", "--version")

    if not docker_available:
        return False, "Docker not installed", None

    # Check if Docker daemon is running
    try:
        result = subprocess.run(["docker", "info"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            return True, docker_version, "Docker daemon running"
        else:
            return False, docker_version, "Docker daemon not running"
    except Exception as e:
        return False, docker_version, f"Cannot check Docker daemon: {e}"


def check_nextflow():
    """Check Nextflow installation and version"""
    nf_available, nf_version = check_command("nextflow", "-version")

    if not nf_available:
        return False, "Nextflow not installed", None

    # Check if version is compatible (>=22.10.0)
    try:
        # Extract version number
        version_line = nf_version.split("\n")[0]
        version_num = version_line.split()[-1]
        major, minor, patch = map(int, version_num.split("."))

        if (major > 22) or (major == 22 and minor >= 10):
            return True, nf_version, "Version compatible"
        else:
            return False, nf_version, f"Version {version_num} too old (need >=22.10.0)"
    except Exception as e:
        return True, nf_version, f"Cannot parse version: {e}"


def verify_pipeline_structure():
    """Verify pipeline directory structure"""
    required_dirs = [
        "workflow",
        "workflow/modules",
        "workflow/conf",
        "containers",
        "docs",
        "scripts",
        "tests/smoke",
        "tests/unit",
        "assets",
        "references",
    ]

    required_files = [
        "workflow/main.nf",
        "workflow/nextflow.config",
        "workflow/conf/base.config",
        "workflow/conf/docker.config",
        "containers/Dockerfile",
        "containers/environment.yml",
        "containers/tool_versions.lock",
        "README.md",
        "LICENSE",
        "CITATION.cff",
    ]

    pipeline_root = Path(__file__).parent.parent

    missing_dirs = []
    missing_files = []

    # Check directories
    for dir_path in required_dirs:
        if not (pipeline_root / dir_path).exists():
            missing_dirs.append(dir_path)

    # Check files
    for file_path in required_files:
        if not (pipeline_root / file_path).exists():
            missing_files.append(file_path)

    return missing_dirs, missing_files


def run_basic_syntax_check():
    """Run basic Nextflow syntax check"""
    pipeline_root = Path(__file__).parent.parent
    workflow_dir = pipeline_root / "workflow"

    try:
        result = subprocess.run(
            ["nextflow", "config", str(workflow_dir)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=workflow_dir,
        )

        if result.returncode == 0:
            return True, "Nextflow config syntax valid"
        else:
            return False, f"Nextflow config error: {result.stderr}"
    except Exception as e:
        return False, f"Cannot run syntax check: {e}"


def check_test_data():
    """Check if test data exists"""
    pipeline_root = Path(__file__).parent.parent

    test_files = [
        "tests/smoke/samplesheet.csv",
        "tests/smoke/tiny_reads_1.fq.gz",
        "tests/smoke/tiny_reads_2.fq.gz",
    ]

    missing_files = []
    for test_file in test_files:
        if not (pipeline_root / test_file).exists():
            missing_files.append(test_file)

    return missing_files


def main():
    parser = argparse.ArgumentParser(description="Verify pipeline installation and dependencies")
    parser.add_argument(
        "--check-docker-pull", action="store_true", help="Also check if Docker image can be pulled"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Show detailed output")

    args = parser.parse_args()

    print("Variant Calling Pipeline - Installation Verification")
    print("=" * 60)

    all_good = True

    # Check Nextflow
    print("\n1. Checking Nextflow...")
    nf_ok, nf_version, nf_status = check_nextflow()
    if nf_ok:
        print(f"   ✓ {nf_version}")
        if nf_status and args.verbose:
            print(f"     {nf_status}")
    else:
        print(f"   ✗ {nf_version}")
        if nf_status:
            print(f"     {nf_status}")
        all_good = False

    # Check Docker
    print("\n2. Checking Docker...")
    docker_ok, docker_version, docker_status = check_docker()
    if docker_ok:
        print(f"   ✓ {docker_version}")
        if docker_status and args.verbose:
            print(f"     {docker_status}")

        # Optionally check Docker image pull
        if args.check_docker_pull:
            print("   Checking Docker image...")
            try:
                result = subprocess.run(
                    ["docker", "pull", "condaforge/mambaforge:4.14.0-0"],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                if result.returncode == 0:
                    print("   ✓ Docker image pull successful")
                else:
                    print("   ⚠ Docker image pull failed (may affect pipeline execution)")
                    if args.verbose:
                        print(f"     {result.stderr}")
            except Exception as e:
                print(f"   ⚠ Cannot test Docker pull: {e}")
    else:
        print(f"   ✗ {docker_version}")
        if docker_status:
            print(f"     {docker_status}")
        all_good = False

    # Check pipeline structure
    print("\n3. Checking pipeline structure...")
    missing_dirs, missing_files = verify_pipeline_structure()

    if not missing_dirs and not missing_files:
        print("   ✓ All required files and directories present")
    else:
        print("   ✗ Missing required components:")
        if missing_dirs:
            print(f"     Missing directories: {', '.join(missing_dirs)}")
        if missing_files:
            print(f"     Missing files: {', '.join(missing_files)}")
        all_good = False

    # Check Nextflow syntax
    print("\n4. Checking Nextflow configuration syntax...")
    syntax_ok, syntax_msg = run_basic_syntax_check()
    if syntax_ok:
        print(f"   ✓ {syntax_msg}")
    else:
        print(f"   ✗ {syntax_msg}")
        all_good = False

    # Check test data
    print("\n5. Checking test data...")
    missing_test_files = check_test_data()
    if not missing_test_files:
        print("   ✓ Test data files present")
    else:
        print("   ⚠ Missing test files (run smoke test setup first):")
        for missing in missing_test_files:
            print(f"     {missing}")
        print("   Run: cd references && bash make_ref_bundle.sh")

    # Summary
    print("\n" + "=" * 60)
    if all_good:
        print("✓ Pipeline installation verification PASSED")
        print("\nYou can now run the pipeline with:")
        print(
            "  nextflow run workflow/main.nf -profile docker --input tests/smoke/samplesheet.csv --outdir results --ref_fasta tests/smoke/ref/test_ref.fa"
        )
    else:
        print("✗ Pipeline installation verification FAILED")
        print("\nPlease fix the issues above before running the pipeline.")
        sys.exit(1)


if __name__ == "__main__":
    main()
