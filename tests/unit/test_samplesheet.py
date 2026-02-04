#!/usr/bin/env python3

"""
Unit tests for samplesheet validation
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path

# Add scripts directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

try:
    from check_samplesheet import check_samplesheet
except ImportError:
    print("Error: Could not import check_samplesheet module")
    sys.exit(1)

class TestSamplesheetValidation(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.temp_dir, "input.csv")
        self.output_file = os.path.join(self.temp_dir, "output.csv")
        
        # Create dummy FASTQ files for testing
        self.read1_file = os.path.join(self.temp_dir, "sample1_R1.fastq.gz")
        self.read2_file = os.path.join(self.temp_dir, "sample1_R2.fastq.gz")
        
        with open(self.read1_file, 'w') as f:
            f.write("@read1\nATCG\n+\nIIII\n")
        with open(self.read2_file, 'w') as f:
            f.write("@read2\nATCG\n+\nIIII\n")
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def create_samplesheet(self, content):
        """Helper to create a samplesheet file"""
        with open(self.input_file, 'w') as f:
            f.write(content)
    
    def test_valid_minimal_samplesheet(self):
        """Test valid samplesheet with minimal required columns"""
        content = f"""sample_id,read1,read2
sample1,{self.read1_file},{self.read2_file}
"""
        self.create_samplesheet(content)
        
        # Should not raise exception
        try:
            check_samplesheet(self.input_file, self.output_file)
            self.assertTrue(os.path.exists(self.output_file))
        except SystemExit:
            self.fail("Valid samplesheet should not raise SystemExit")
    
    def test_valid_full_samplesheet(self):
        """Test valid samplesheet with all columns"""
        content = f"""sample_id,read1,read2,platform,library,lane
sample1,{self.read1_file},{self.read2_file},ILLUMINA,LIB1,1
sample2,{self.read1_file},{self.read2_file},ILLUMINA,LIB2,2
"""
        self.create_samplesheet(content)
        
        try:
            check_samplesheet(self.input_file, self.output_file)
            self.assertTrue(os.path.exists(self.output_file))
        except SystemExit:
            self.fail("Valid samplesheet should not raise SystemExit")
    
    def test_missing_header(self):
        """Test samplesheet without header"""
        content = f"""sample1,{self.read1_file},{self.read2_file}
"""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_missing_required_column(self):
        """Test samplesheet missing required column"""
        content = f"""sample_id,read1
sample1,{self.read1_file}
"""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_empty_sample_id(self):
        """Test samplesheet with empty sample_id"""
        content = f""",read1,read2
,{self.read1_file},{self.read2_file}
"""
        # Add header line
        content = "sample_id,read1,read2\n" + content
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_duplicate_sample_id(self):
        """Test samplesheet with duplicate sample IDs"""
        content = f"""sample_id,read1,read2
sample1,{self.read1_file},{self.read2_file}
sample1,{self.read1_file},{self.read2_file}
"""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_missing_fastq_file(self):
        """Test samplesheet with non-existent FASTQ file"""
        content = """sample_id,read1,read2
sample1,/nonexistent/file_R1.fastq.gz,/nonexistent/file_R2.fastq.gz
"""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_invalid_file_extension(self):
        """Test samplesheet with invalid file extensions"""
        # Create a file with wrong extension
        wrong_ext_file = os.path.join(self.temp_dir, "sample.txt")
        with open(wrong_ext_file, 'w') as f:
            f.write("test")
        
        content = f"""sample_id,read1,read2
sample1,{wrong_ext_file},{self.read2_file}
"""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_valid_file_extensions(self):
        """Test all valid file extensions"""
        valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        
        for i, ext in enumerate(valid_extensions):
            read1_file = os.path.join(self.temp_dir, f"sample{i}_R1{ext}")
            read2_file = os.path.join(self.temp_dir, f"sample{i}_R2{ext}")
            
            with open(read1_file, 'w') as f:
                f.write("test")
            with open(read2_file, 'w') as f:
                f.write("test")
            
            content = f"""sample_id,read1,read2
sample{i},{read1_file},{read2_file}
"""
            input_file = os.path.join(self.temp_dir, f"input{i}.csv")
            output_file = os.path.join(self.temp_dir, f"output{i}.csv")
            
            with open(input_file, 'w') as f:
                f.write(content)
            
            try:
                check_samplesheet(input_file, output_file)
                self.assertTrue(os.path.exists(output_file))
            except SystemExit:
                self.fail(f"Valid extension {ext} should not raise SystemExit")
    
    def test_default_values(self):
        """Test that default values are properly assigned"""
        content = f"""sample_id,read1,read2
sample1,{self.read1_file},{self.read2_file}
"""
        self.create_samplesheet(content)
        
        try:
            check_samplesheet(self.input_file, self.output_file)
            
            # Read the output file and check defaults
            with open(self.output_file, 'r') as f:
                lines = f.readlines()
                
            # Should have header + 1 data line
            self.assertEqual(len(lines), 2)
            
            # Check header
            header = lines[0].strip().split(',')
            expected_header = ["sample_id", "read1", "read2", "platform", "library", "lane"]
            self.assertEqual(header, expected_header)
            
            # Check data line has defaults
            data = lines[1].strip().split(',')
            self.assertEqual(data[3], "ILLUMINA")  # Default platform
            self.assertEqual(data[4], "sample1")   # Default library (sample_id)
            self.assertEqual(data[5], "1")         # Default lane
            
        except SystemExit:
            self.fail("Valid samplesheet should not raise SystemExit")
    
    def test_whitespace_handling(self):
        """Test that whitespace is properly stripped"""
        content = f"""sample_id,read1,read2,platform,library,lane
 sample1 , {self.read1_file} , {self.read2_file} , ILLUMINA , LIB1 , 1 
"""
        self.create_samplesheet(content)
        
        try:
            check_samplesheet(self.input_file, self.output_file)
            self.assertTrue(os.path.exists(self.output_file))
            
            # Read output and verify whitespace was stripped
            with open(self.output_file, 'r') as f:
                lines = f.readlines()
            
            data = lines[1].strip().split(',')
            self.assertEqual(data[0], "sample1")  # No leading/trailing spaces
            
        except SystemExit:
            self.fail("Samplesheet with whitespace should not raise SystemExit")
    
    def test_empty_file(self):
        """Test completely empty samplesheet"""
        content = ""
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)
    
    def test_header_only(self):
        """Test samplesheet with header but no data"""
        content = "sample_id,read1,read2\n"
        self.create_samplesheet(content)
        
        with self.assertRaises(SystemExit):
            check_samplesheet(self.input_file, self.output_file)

class TestSamplesheetParsing(unittest.TestCase):
    """Additional tests for edge cases and parsing logic"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_case_insensitive_headers(self):
        """Test that headers are case insensitive"""
        # Create dummy files
        read1_file = os.path.join(self.temp_dir, "R1.fastq.gz")
        read2_file = os.path.join(self.temp_dir, "R2.fastq.gz")
        with open(read1_file, 'w') as f:
            f.write("test")
        with open(read2_file, 'w') as f:
            f.write("test")
        
        content = f"""SAMPLE_ID,READ1,READ2
sample1,{read1_file},{read2_file}
"""
        input_file = os.path.join(self.temp_dir, "input.csv")
        output_file = os.path.join(self.temp_dir, "output.csv")
        
        with open(input_file, 'w') as f:
            f.write(content)
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(os.path.exists(output_file))
        except SystemExit:
            self.fail("Case insensitive headers should be accepted")

def run_tests():
    """Run all tests"""
    unittest.main(verbosity=2)

if __name__ == '__main__':
    run_tests()