#!/usr/bin/env python3

"""
Unit tests for VCF quality control script
"""

import unittest
import tempfile
import os
import sys
import gzip
import json
from pathlib import Path

# Add scripts directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

try:
    from vcf_qc import VCFQualityControl
except ImportError:
    print("Error: Could not import vcf_qc module")
    sys.exit(1)

class TestVCFQualityControl(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.vcf_file = os.path.join(self.temp_dir, "test.vcf")
        self.vcf_gz_file = os.path.join(self.temp_dir, "test.vcf.gz")
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def create_test_vcf(self, compressed=False):
        """Create a test VCF file with known variants"""
        vcf_content = """##fileformat=VCFv4.2
##reference=test_ref.fa
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by Depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="FisherStrand">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping Quality">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand Odds Ratio">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	rs1	A	T	30.0	PASS	DP=20;QD=1.5;FS=0.5;MQ=40.0;SOR=1.2	GT:GQ:DP	0/1:30:20
chr1	200	rs2	T	C	50.0	PASS	DP=25;QD=2.0;FS=1.0;MQ=45.0;SOR=0.8	GT:GQ:DP	1/1:40:25
chr1	300	rs3	C	G	20.0	LowQual	DP=15;QD=1.0;FS=2.0;MQ=30.0;SOR=2.5	GT:GQ:DP	0/1:20:15
chr1	400	rs4	G	A	60.0	PASS	DP=30;QD=2.5;FS=0.2;MQ=50.0;SOR=0.5	GT:GQ:DP	0/1:50:30
chr1	500	rs5	AT	A	40.0	PASS	DP=22;QD=1.8;FS=1.5;MQ=42.0;SOR=1.0	GT:GQ:DP	0/1:35:22
"""
        
        if compressed:
            with gzip.open(self.vcf_gz_file, 'wt') as f:
                f.write(vcf_content)
            return self.vcf_gz_file
        else:
            with open(self.vcf_file, 'w') as f:
                f.write(vcf_content)
            return self.vcf_file
    
    def create_minimal_vcf(self):
        """Create minimal valid VCF with header only"""
        vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""
        with open(self.vcf_file, 'w') as f:
            f.write(vcf_content)
        return self.vcf_file
    
    def create_invalid_vcf(self):
        """Create invalid VCF file"""
        vcf_content = """This is not a VCF file
chr1	100	rs1	A	T	30.0
"""
        with open(self.vcf_file, 'w') as f:
            f.write(vcf_content)
        return self.vcf_file
    
    def test_vcf_validation_valid_file(self):
        """Test VCF validation with valid file"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        issues = qc.validate_vcf()
        self.assertEqual(len(issues), 0, "Valid VCF should have no validation issues")
    
    def test_vcf_validation_missing_file(self):
        """Test VCF validation with missing file"""
        missing_file = os.path.join(self.temp_dir, "missing.vcf")
        qc = VCFQualityControl(missing_file, self.temp_dir)
        
        issues = qc.validate_vcf()
        self.assertGreater(len(issues), 0, "Missing VCF should have validation issues")
        self.assertIn("does not exist", issues[0])
    
    def test_vcf_validation_invalid_format(self):
        """Test VCF validation with invalid format"""
        vcf_file = self.create_invalid_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        issues = qc.validate_vcf()
        self.assertGreater(len(issues), 0, "Invalid VCF should have validation issues")
    
    def test_vcf_validation_empty_file(self):
        """Test VCF validation with empty file"""
        with open(self.vcf_file, 'w') as f:
            pass  # Create empty file
        
        qc = VCFQualityControl(self.vcf_file, self.temp_dir)
        issues = qc.validate_vcf()
        self.assertGreater(len(issues), 0, "Empty VCF should have validation issues")
        self.assertIn("empty", issues[0])
    
    def test_compressed_vcf_reading(self):
        """Test reading compressed VCF files"""
        vcf_file = self.create_test_vcf(compressed=True)
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        # Should be able to analyze compressed VCF
        qc.analyze_vcf()
        
        # Check that variants were found
        self.assertGreater(qc.stats['total_variants'], 0)
        self.assertEqual(qc.stats['total_variants'], 5)  # 5 variants in test data
    
    def test_variant_type_classification(self):
        """Test variant type classification"""
        qc = VCFQualityControl(self.vcf_file, self.temp_dir)
        
        # Test SNP transitions
        self.assertEqual(qc.classify_variant("A", "G"), "transition")
        self.assertEqual(qc.classify_variant("C", "T"), "transition")
        self.assertEqual(qc.classify_variant("G", "A"), "transition")
        self.assertEqual(qc.classify_variant("T", "C"), "transition")
        
        # Test SNP transversions
        self.assertEqual(qc.classify_variant("A", "C"), "transversion")
        self.assertEqual(qc.classify_variant("A", "T"), "transversion")
        self.assertEqual(qc.classify_variant("G", "C"), "transversion")
        self.assertEqual(qc.classify_variant("G", "T"), "transversion")
        
        # Test indels
        self.assertEqual(qc.classify_variant("AT", "A"), "indel")  # Deletion
        self.assertEqual(qc.classify_variant("A", "AT"), "indel")  # Insertion
        
        # Test complex variants
        self.assertEqual(qc.classify_variant("ATG", "CTA"), "complex")
    
    def test_info_field_parsing(self):
        """Test INFO field parsing"""
        qc = VCFQualityControl(self.vcf_file, self.temp_dir)
        
        # Test typical INFO field
        info_str = "DP=20;QD=1.5;FS=0.5;MQ=40.0;PASS"
        info_dict = qc.parse_info_field(info_str)
        
        self.assertEqual(info_dict['DP'], 20)
        self.assertEqual(info_dict['QD'], 1.5)
        self.assertEqual(info_dict['FS'], 0.5)
        self.assertEqual(info_dict['MQ'], 40.0)
        self.assertTrue(info_dict['PASS'])
        
        # Test empty INFO field
        empty_info = qc.parse_info_field(".")
        self.assertEqual(len(empty_info), 0)
    
    def test_vcf_analysis_statistics(self):
        """Test VCF analysis and statistics calculation"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        qc.analyze_vcf()
        
        # Check basic counts
        self.assertEqual(qc.stats['total_variants'], 5)
        self.assertEqual(qc.stats['snps'], 4)  # 4 SNPs in test data
        self.assertEqual(qc.stats['indels'], 1)  # 1 indel in test data
        self.assertEqual(qc.stats['passed'], 4)  # 4 variants passed filters
        self.assertEqual(qc.stats['filtered'], 1)  # 1 variant filtered
        
        # Check variant type classification
        self.assertEqual(qc.stats['transitions'], 2)  # A>T and T>C transitions
        self.assertEqual(qc.stats['transversions'], 2)  # C>G and G>A transversions
        
        # Check Ti/Tv ratio calculation
        expected_ti_tv = 2.0 / 2.0  # 2 transitions / 2 transversions = 1.0
        self.assertAlmostEqual(qc.stats['ti_tv_ratio'], expected_ti_tv, places=2)
        
        # Check that quality data was collected
        self.assertGreater(len(qc.quality_data['QUAL']), 0)
        self.assertGreater(len(qc.quality_data['DP']), 0)
    
    def test_json_output(self):
        """Test JSON statistics output"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        qc.analyze_vcf()
        json_stats = qc.save_json_stats()
        
        # Check that JSON output contains expected fields
        self.assertIn('total_variants', json_stats)
        self.assertIn('snps', json_stats)
        self.assertIn('indels', json_stats)
        self.assertIn('ti_tv_ratio', json_stats)
        self.assertIn('quality_metrics', json_stats)
        self.assertIn('variant_types', json_stats)
        
        # Verify JSON file was created
        json_file = Path(self.temp_dir) / "test_qc_stats.json"
        self.assertTrue(json_file.exists())
        
        # Verify JSON file can be loaded
        with open(json_file, 'r') as f:
            loaded_stats = json.load(f)
        self.assertEqual(loaded_stats['total_variants'], json_stats['total_variants'])
    
    def test_report_generation(self):
        """Test text report generation"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        qc.analyze_vcf()
        qc.generate_report()
        
        # Check that report file was created
        report_file = Path(self.temp_dir) / "test_qc_report.txt"
        self.assertTrue(report_file.exists())
        
        # Check report content
        with open(report_file, 'r') as f:
            report_content = f.read()
        
        self.assertIn("VCF Quality Control Report", report_content)
        self.assertIn("SUMMARY STATISTICS", report_content)
        self.assertIn("QUALITY METRICS", report_content)
        self.assertIn("Total variants:", report_content)
        self.assertIn("Ti/Tv ratio:", report_content)
    
    def test_minimal_vcf_handling(self):
        """Test handling of minimal VCF with header only"""
        vcf_file = self.create_minimal_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        # Validation should pass (has proper header)
        issues = qc.validate_vcf()
        self.assertEqual(len(issues), 0)
        
        # Analysis should handle zero variants gracefully
        qc.analyze_vcf()
        self.assertEqual(qc.stats['total_variants'], 0)
        self.assertEqual(qc.stats['ti_tv_ratio'], 0.0)
    
    def test_filter_counting(self):
        """Test filter status counting"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        qc.analyze_vcf()
        
        # Check filter counts
        self.assertIn('LowQual', qc.filter_counts)
        self.assertEqual(qc.filter_counts['LowQual'], 1)
    
    def test_quality_metrics_calculation(self):
        """Test quality metrics calculation"""
        vcf_file = self.create_test_vcf()
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        qc.analyze_vcf()
        
        # Check that quality metrics were calculated
        self.assertGreater(qc.stats['quality_metrics']['mean_qual'], 0)
        self.assertGreater(qc.stats['quality_metrics']['median_qual'], 0)
        self.assertGreater(qc.stats['quality_metrics']['mean_dp'], 0)
        self.assertGreater(qc.stats['quality_metrics']['median_dp'], 0)
        
        # Verify specific values from test data
        expected_mean_qual = (30.0 + 50.0 + 20.0 + 60.0 + 40.0) / 5  # 40.0
        self.assertAlmostEqual(qc.stats['quality_metrics']['mean_qual'], expected_mean_qual, places=1)

class TestVCFQCIntegration(unittest.TestCase):
    """Integration tests for VCF QC functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_ti_tv_ratio_validation(self):
        """Test Ti/Tv ratio calculation for different variant sets"""
        # Create VCF with known Ti/Tv ratio
        vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr1	200	.	C	T	30	PASS	.
chr1	300	.	G	A	30	PASS	.
chr1	400	.	T	C	30	PASS	.
chr1	500	.	A	C	30	PASS	.
chr1	600	.	A	T	30	PASS	.
"""
        
        vcf_file = os.path.join(self.temp_dir, "test_titv.vcf")
        with open(vcf_file, 'w') as f:
            f.write(vcf_content)
        
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        qc.analyze_vcf()
        
        # Should have 4 transitions and 2 transversions
        self.assertEqual(qc.stats['transitions'], 4)
        self.assertEqual(qc.stats['transversions'], 2)
        self.assertAlmostEqual(qc.stats['ti_tv_ratio'], 2.0, places=2)
    
    def test_edge_cases(self):
        """Test edge cases and error handling"""
        # Test with malformed VCF line
        vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr1	invalid_line_with_few_fields
chr1	300	.	C	T	40	PASS	.
"""
        
        vcf_file = os.path.join(self.temp_dir, "test_malformed.vcf")
        with open(vcf_file, 'w') as f:
            f.write(vcf_content)
        
        qc = VCFQualityControl(vcf_file, self.temp_dir)
        
        # Should handle malformed lines gracefully
        qc.analyze_vcf()
        self.assertEqual(qc.stats['total_variants'], 2)  # Should skip malformed line

def run_tests():
    """Run all tests"""
    unittest.main(verbosity=2)

if __name__ == '__main__':
    run_tests()