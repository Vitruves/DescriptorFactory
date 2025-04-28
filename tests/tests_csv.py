import pytest
import os
import csv
import subprocess
import tempfile
import shutil
import pandas as pd
from pathlib import Path


@pytest.fixture
def desfact_executable():
    return "./build/desfact"


@pytest.fixture
def temp_dir():
    dir_path = tempfile.mkdtemp()
    yield dir_path
    shutil.rmtree(dir_path)


@pytest.fixture
def valid_smiles_csv(temp_dir):
    csv_path = Path(temp_dir) / "valid_smiles.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Property"])
        writer.writerow(["1", "CCO", "Ethanol"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid"])
        writer.writerow(["3", "c1ccccc1", "Benzene"])
        writer.writerow(["4", "C1CCCCC1", "Cyclohexane"])
        writer.writerow(["5", "CC(C)C", "Isobutane"])
    
    return csv_path


@pytest.fixture
def custom_delimiter_csv(temp_dir):
    csv_path = Path(temp_dir) / "custom_delimiter.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(["ID", "SMILES", "Property"])
        writer.writerow(["1", "CCO", "Ethanol"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid"])
    
    return csv_path


@pytest.fixture
def quoted_fields_csv(temp_dir):
    csv_path = Path(temp_dir) / "quoted_fields.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Name with, comma"])
        writer.writerow(["1", "CCO", "Ethanol, pure"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid, glacial"])
    
    return csv_path


@pytest.fixture
def newlines_in_fields_csv(temp_dir):
    csv_path = Path(temp_dir) / "newlines_in_fields.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Multi-line\nDescription"])
        writer.writerow(["1", "CCO", "Ethanol\nAlso called drinking alcohol"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid\nFound in vinegar"])
    
    return csv_path


@pytest.fixture
def missing_smiles_csv(temp_dir):
    csv_path = Path(temp_dir) / "missing_smiles.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "Structure", "Property"])
        writer.writerow(["1", "CCO", "Ethanol"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid"])
    
    return csv_path


@pytest.fixture
def invalid_smiles_csv(temp_dir):
    csv_path = Path(temp_dir) / "invalid_smiles.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Property"])
        writer.writerow(["1", "CCO", "Ethanol"])
        writer.writerow(["2", "INVALID", "Invalid compound"])
        writer.writerow(["3", "c1ccccc1", "Benzene"])
    
    return csv_path


@pytest.fixture
def no_header_csv(temp_dir):
    csv_path = Path(temp_dir) / "no_header.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["1", "CCO", "Ethanol"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid"])
    
    return csv_path


@pytest.fixture
def empty_fields_csv(temp_dir):
    csv_path = Path(temp_dir) / "empty_fields.csv"
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Property"])
        writer.writerow(["1", "CCO", ""])
        writer.writerow(["2", "", "Acetic acid"])
        writer.writerow(["3", "c1ccccc1", ""])
    
    return csv_path


@pytest.fixture
def malformed_csv(temp_dir):
    csv_path = Path(temp_dir) / "malformed.csv"
    
    with open(csv_path, 'w', newline='') as f:
        f.write("ID,SMILES,Property\n")
        f.write("1,CCO,Ethanol\n")
        f.write("2,CC(=O)O,Acetic acid\n")
        f.write("3,c1ccccc1\n")  # Missing field
        f.write("4,C1CCCCC1,Cyclohexane\n")
    
    return csv_path


def test_all_descriptors(desfact_executable, valid_smiles_csv, temp_dir):
    output_path = Path(temp_dir) / "output.csv"
    
    cmd = [
        desfact_executable,
        "-i", str(valid_smiles_csv),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    df_in = pd.read_csv(valid_smiles_csv)
    df_out = pd.read_csv(output_path)
    
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    assert df_out.shape[1] > df_in.shape[1]  # Should have additional descriptor columns


def test_preservation_of_order(desfact_executable, valid_smiles_csv, temp_dir):
    output_path = Path(temp_dir) / "output.csv"
    
    cmd = [
        desfact_executable,
        "-i", str(valid_smiles_csv),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    df_in = pd.read_csv(valid_smiles_csv)
    df_out = pd.read_csv(output_path)
    
    # Check that the original rows are preserved in the same order
    for i in range(len(df_in)):
        for col in df_in.columns:
            if col != "SMILES":  # Skip SMILES column which might be canonicalized
                assert df_in.loc[i, col] == df_out.loc[i, col]


def test_large_file(desfact_executable, temp_dir):
    large_csv_path = Path(temp_dir) / "large.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create a large CSV with 1000 rows
    with open(large_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Property"])
        
        for i in range(1000):
            # Alternate between a few valid SMILES
            if i % 5 == 0:
                smiles = "CCO"
            elif i % 5 == 1:
                smiles = "CC(=O)O"
            elif i % 5 == 2:
                smiles = "c1ccccc1"
            elif i % 5 == 3:
                smiles = "C1CCCCC1"
            else:
                smiles = "CC(C)C"
                
            writer.writerow([str(i+1), smiles, f"Compound {i+1}"])
    
    cmd = [
        desfact_executable,
        "-i", str(large_csv_path),
        "-o", str(output_path),
        "-d", "all",
        "-b", "100",  # Process in batches of 100
        "-t", "4"     # Use 4 threads
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    df_in = pd.read_csv(large_csv_path)
    df_out = pd.read_csv(output_path)
    
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    assert df_out.shape[1] > df_in.shape[1]


def test_csv_integrity_with_all_descriptors(desfact_executable, valid_smiles_csv, temp_dir):
    output_path = Path(temp_dir) / "output.csv"
    
    cmd = [
        desfact_executable,
        "-i", str(valid_smiles_csv),
        "-o", str(output_path),
        "-d", "all",
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    assert os.path.exists(output_path)
    
    df_in = pd.read_csv(valid_smiles_csv)
    df_out = pd.read_csv(output_path)
    
    # Test basic integrity
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    
    # Test preservation of non-SMILES values
    for i in range(len(df_in)):
        for col in df_in.columns:
            if col != "SMILES":  # SMILES might be canonicalized
                assert df_in.loc[i, col] == df_out.loc[i, col]
    
    # Test no missing values in descriptor columns
    descriptor_cols = [col for col in df_out.columns if col not in df_in.columns]
    assert len(descriptor_cols) > 0
    for col in descriptor_cols:
        assert not df_out[col].isna().any(), f"Found missing values in descriptor column: {col}"


def test_csv_encoding_with_all_descriptors(desfact_executable, temp_dir):
    # Create CSV with special characters
    csv_path = Path(temp_dir) / "special_chars.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Name"])
        writer.writerow(["1", "CCO", "Ethanol-α"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid-β"])
        writer.writerow(["3", "c1ccccc1", "Benzène"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check encoding is preserved
    df_in = pd.read_csv(csv_path, encoding='utf-8')
    df_out = pd.read_csv(output_path, encoding='utf-8')
    
    for i in range(len(df_in)):
        assert df_in.loc[i, "Name"] == df_out.loc[i, "Name"]


def test_get_all(desfact_executable):
    """Get a valid descriptor name to use in tests"""
    cmd = [desfact_executable, "-l"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Parse the output to find available descriptors
    output_lines = result.stdout.strip().split('\n')
    descriptor_lines = [line.strip() for line in output_lines if ':' in line]
    
    # Extract the descriptor names
    valid_descriptors = []
    for line in descriptor_lines:
        parts = line.split(':')
        if len(parts) >= 2:
            descriptor_name = parts[0].strip()
            valid_descriptors.append(descriptor_name)
    
    assert len(valid_descriptors) > 0, "No valid descriptors found"
    
    # Use the first available descriptor
    return valid_descriptors[0]


@pytest.fixture
def valid_descriptor(desfact_executable):
    """Fixture to provide a valid descriptor name"""
    return test_get_all(desfact_executable)


def test_complex_csv_structure_with_all_descriptors(desfact_executable, temp_dir):
    """Test handling of complex CSV structure with multiple data types."""
    csv_path = Path(temp_dir) / "complex.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create a CSV with various data types
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["ID", "SMILES", "Name", "Numbers", "Boolean", "Date"])
        writer.writerow(["1", "CCO", "Ethanol", "123.45", "TRUE", "2023-01-01"])
        writer.writerow(["2", "CC(=O)O", "Acetic acid", "-456.78", "FALSE", "2023-02-15"])
        writer.writerow(["3", "c1ccccc1", "Benzene", "0", "TRUE", "2023-03-30"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all",
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    # Verify data integrity for all original columns
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    
    # Check each non-SMILES column is preserved exactly
    for col in df_in.columns:
        if col != "SMILES":  # SMILES might be canonicalized
            for i in range(len(df_in)):
                assert str(df_in.loc[i, col]) == str(df_out.loc[i, col])
    
    # Verify descriptor columns were added
    assert df_out.shape[1] > df_in.shape[1]


def test_escaped_csv_with_all_descriptors(desfact_executable, temp_dir):
    """Test CSV with escaped delimiters and quotes."""
    csv_path = Path(temp_dir) / "escaped.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create a CSV with escaped characters
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Description"])
        writer.writerow(["1", "CCO", "Contains, a comma"])
        writer.writerow(["2", "CC(=O)O", "Has \"quotes\" inside"])
        writer.writerow(["3", "c1ccccc1", "Multiple lines\nIn one cell"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Need to use the same quoting style when reading
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    # Verify all special characters are preserved
    for i in range(len(df_in)):
        if "comma" in df_in.loc[i, "Description"]:
            assert "," in df_out.loc[i, "Description"]
        if "quotes" in df_in.loc[i, "Description"]:
            assert "\"" in df_out.loc[i, "Description"]
        if "Multiple lines" in df_in.loc[i, "Description"]:
            assert "\n" in df_out.loc[i, "Description"]


def test_multiline_csv_records_with_all_descriptors(desfact_executable, temp_dir):
    """Test CSV with multiline records."""
    csv_path = Path(temp_dir) / "multiline.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "LongDescription"])
        writer.writerow(["1", "CCO", "Line 1\nLine 2\nLine 3"])
        writer.writerow(["2", "CC(=O)O", "Another\nmultiline\ndescription"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check if output exists and has at least as many descriptors as expected
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    
    # Verify newlines in descriptions are preserved
    for i in range(len(df_in)):
        desc_in = df_in.loc[i, "LongDescription"]
        desc_out = df_out.loc[i, "LongDescription"]
        assert desc_in.count('\n') == desc_out.count('\n')


def test_csv_handling_with_mixed_numeric_types(desfact_executable, temp_dir):
    """Test handling of CSVs containing mixed numeric types (integers, floats, scientific notation)."""
    csv_path = Path(temp_dir) / "mixed_numeric.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with mixed numeric types
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Integer", "Float", "Scientific", "Negative"])
        writer.writerow(["1", "CCO", "123", "45.67", "1.23e-4", "-789"])
        writer.writerow(["2", "CC(=O)O", "0", "0.0001", "6.022e23", "-0.5"])
        writer.writerow(["3", "c1ccccc1", "9999", "3.14159", "1.602e-19", "-99.99"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Verify numeric values are preserved precisely
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    numeric_cols = ["Integer", "Float", "Scientific", "Negative"]
    for col in numeric_cols:
        for i in range(len(df_in)):
            # Convert to string to ensure format preservation
            assert str(df_in.loc[i, col]) == str(df_out.loc[i, col])


def test_csv_handling_special_characters(desfact_executable, temp_dir):
    """Test handling of CSVs with special characters including Unicode."""
    csv_path = Path(temp_dir) / "special_chars.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with special characters
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "SpecialChars"])
        writer.writerow(["1", "CCO", "!@#$%^&*()_+{}[]|\\:;\"'<>,.?/"])
        writer.writerow(["2", "CC(=O)O", "±§½¾¿¡œæßðøþ"])
        writer.writerow(["3", "c1ccccc1", "你好こんにちは안녕하세요"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Verify special characters are preserved
    df_in = pd.read_csv(csv_path, encoding='utf-8', quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, encoding='utf-8', quoting=csv.QUOTE_ALL)
    
    for i in range(len(df_in)):
        assert df_in.loc[i, "SpecialChars"] == df_out.loc[i, "SpecialChars"]


def test_csv_quoting_styles(desfact_executable, temp_dir):
    """Test handling of different CSV quoting styles."""
    # Test different quoting styles
    for quoting_style in [csv.QUOTE_ALL, csv.QUOTE_MINIMAL, csv.QUOTE_NONNUMERIC]:
        csv_path = Path(temp_dir) / f"quoting_{quoting_style}.csv"
        output_path = Path(temp_dir) / f"output_{quoting_style}.csv"
        
        # Create CSV with specific quoting style
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f, quoting=quoting_style)
            writer.writerow(["ID", "SMILES", "Text", "Number"])
            writer.writerow(["1", "CCO", "With, comma", "123"])
            writer.writerow(["2", "CC(=O)O", "With \"quotes\"", "456"])
        
        cmd = [
            desfact_executable,
            "-i", str(csv_path),
            "-o", str(output_path),
            "-d", "all"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        
        # Verify content is preserved
        df_in = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_path)
        
        for i in range(len(df_in)):
            assert df_in.loc[i, "Text"] == df_out.loc[i, "Text"]
            assert str(df_in.loc[i, "Number"]) == str(df_out.loc[i, "Number"])


def test_csv_custom_escape_character(desfact_executable, temp_dir):
    """Test handling of CSV with custom escape characters."""
    csv_path = Path(temp_dir) / "escaped.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Manually create CSV with backslash escapes
    with open(csv_path, 'w', newline='') as f:
        f.write('ID,SMILES,Text\n')
        f.write('1,CCO,"Text with \\"escaped quotes\\" inside"\n')
        f.write('2,CC(=O)O,"Line\\nwith\\nescaped\\nnewlines"\n')
        f.write('3,c1ccccc1,"Backslash\\\\character"\n')
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check if output exists and escaped sequences are preserved
    df_in = pd.read_csv(csv_path, escapechar='\\')
    df_out = pd.read_csv(output_path, escapechar='\\')
    
    # Verify escaped characters are handled correctly
    for i in range(len(df_in)):
        if i == 0:
            assert '"escaped quotes"' in df_out.loc[i, "Text"]
        elif i == 1:
            assert df_out.loc[i, "Text"].count('\\n') == 3
        elif i == 2:
            assert '\\\\' in df_out.loc[i, "Text"]


def test_csv_header_preservation(desfact_executable, temp_dir):
    """Test preservation of CSV headers including special characters and whitespace."""
    csv_path = Path(temp_dir) / "headers.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with unusual headers
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Column With Spaces", "Column-With-Hyphens", "Column.With.Dots", "Column_With_Underscores"])
        writer.writerow(["1", "CCO", "a", "b", "c", "d"])
        writer.writerow(["2", "CC(=O)O", "e", "f", "g", "h"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check header preservation
    with open(csv_path, 'r') as f:
        in_headers = next(csv.reader(f))
    
    with open(output_path, 'r') as f:
        out_headers = next(csv.reader(f))
    
    # Verify all original headers exist in output headers
    for header in in_headers:
        assert header in out_headers


def test_csv_empty_file_handling(desfact_executable, temp_dir):
    """Test handling of empty CSV files."""
    # Test various empty file scenarios
    scenarios = [
        ("empty.csv", ""),  # Completely empty file
        ("header_only.csv", "ID,SMILES,Property\n"),  # Header row only
        ("no_data.csv", "ID,SMILES,Property\n\n\n")  # Header row with blank lines
    ]
    
    for filename, content in scenarios:
        csv_path = Path(temp_dir) / filename
        output_path = Path(temp_dir) / f"output_{filename}"
        
        # Create the test file
        with open(csv_path, 'w') as f:
            f.write(content)
        
        cmd = [
            desfact_executable,
            "-i", str(csv_path),
            "-o", str(output_path),
            "-d", "all"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Program should either handle empty files gracefully or exit with informative error
        if result.returncode == 0:
            assert os.path.exists(output_path), "Output file should be created"
        else:
            assert "empty" in result.stderr.lower() or "no data" in result.stderr.lower() or "error" in result.stderr.lower()


def test_csv_binary_data_fields(desfact_executable, temp_dir):
    """Test handling of CSV with binary data encoded as base64."""
    csv_path = Path(temp_dir) / "binary_data.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Sample base64 encoded binary data
    base64_data = "SGVsbG8gV29ybGQh"  # "Hello World!" in base64
    base64_image = "R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7"  # 1x1 transparent GIF
    
    # Create CSV with base64 encoded fields
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "BinaryData"])
        writer.writerow(["1", "CCO", base64_data])
        writer.writerow(["2", "CC(=O)O", base64_image])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Verify base64 strings are preserved exactly
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    assert df_in.loc[0, "BinaryData"] == df_out.loc[0, "BinaryData"]
    assert df_in.loc[1, "BinaryData"] == df_out.loc[1, "BinaryData"]


def test_csv_very_wide_file(desfact_executable, temp_dir):
    """Test handling of CSV files with a large number of columns."""
    csv_path = Path(temp_dir) / "wide.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create a CSV with 100 columns
    header = ["ID", "SMILES"] + [f"Col{i}" for i in range(98)]
    data_row1 = ["1", "CCO"] + [f"value1_{i}" for i in range(98)]
    data_row2 = ["2", "CC(=O)O"] + [f"value2_{i}" for i in range(98)]
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data_row1)
        writer.writerow(data_row2)
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Verify all columns are preserved
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    assert len(df_in.columns) == 100
    assert all(col in df_out.columns for col in df_in.columns)
    
    # Check data integrity for all original columns
    for col in df_in.columns:
        if col != "SMILES":  # SMILES might be canonicalized
            for i in range(len(df_in)):
                assert df_in.loc[i, col] == df_out.loc[i, col]


def test_error_messages(desfact_executable, temp_dir):
    """Test that appropriate error messages are shown for various error conditions."""
    
    # Test nonexistent file
    nonexistent_path = Path(temp_dir) / "nonexistent.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    cmd = [
        desfact_executable,
        "-i", str(nonexistent_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    error_output = result.stderr.lower() + result.stdout.lower()
    assert any(term in error_output for term in ["no such file", "cannot open", "not found", "error"])
    
    # Test missing SMILES column
    no_smiles_path = Path(temp_dir) / "no_smiles.csv"
    with open(no_smiles_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "Structure", "Property"])
        writer.writerow(["1", "CCO", "Ethanol"])
    
    cmd = [
        desfact_executable,
        "-i", str(no_smiles_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    error_output = result.stderr.lower() + result.stdout.lower()
    assert any(term in error_output for term in ["smiles", "column", "not found", "error"])


def test_csv_output_integrity_with_invalid_smiles(desfact_executable, temp_dir):
    """Test output integrity when some SMILES are invalid."""
    csv_path = Path(temp_dir) / "mixed_valid_invalid.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with mix of valid and invalid SMILES
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Name"])
        writer.writerow(["1", "CCO", "Valid 1"])
        writer.writerow(["2", "XXXXX", "Invalid"])  # Invalid SMILES
        writer.writerow(["3", "c1ccccc1", "Valid 2"])
        writer.writerow(["4", "CC(=O)O", "Valid 3"])
        writer.writerow(["5", "C1C(C)C1(CC)C", "Invalid Structure"])  # Invalid structure
        writer.writerow(["6", "c1ccccc1O", "Valid 4"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all",
        "--verbose"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check if it continues with warnings or exits with error
    if result.returncode == 0:
        # If it continues, check that valid rows are processed correctly
        df_in = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_path)
        
        # All rows should be present
        assert len(df_in) == len(df_out)
        
        # Check that original columns are preserved
        assert all(col in df_out.columns for col in df_in.columns)
        
        # Check for descriptor columns
        assert df_out.shape[1] > df_in.shape[1]
        
        # Non-SMILES data for valid molecules should be preserved
        valid_indices = [0, 2, 3, 5]  # Indices of rows with valid SMILES
        for i in valid_indices:
            assert df_in.loc[i, "Name"] == df_out.loc[i, "Name"]
            
        # Check if warning/error messages were generated for invalid SMILES
        assert "invalid" in result.stderr.lower() or "invalid" in result.stdout.lower()
    else:
        # If it exits with error, check for appropriate error message
        assert "invalid" in result.stderr.lower() or "invalid" in result.stdout.lower()


def test_csv_row_integrity_with_missing_values(desfact_executable, temp_dir):
    """Test that row integrity is maintained with missing values."""
    csv_path = Path(temp_dir) / "missing_values.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with missing values
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Name", "Property1", "Property2"])
        writer.writerow(["1", "CCO", "Ethanol", "1.0", ""])  # Missing Property2
        writer.writerow(["2", "CC(=O)O", "", "2.0", "b"])    # Missing Name
        writer.writerow(["3", "c1ccccc1", "Benzene", "", ""]) # Missing both properties
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    # Check row count
    assert len(df_in) == len(df_out)
    
    # Check empty values are preserved
    for i in range(len(df_in)):
        for col in ["Name", "Property1", "Property2"]:
            if pd.isna(df_in.loc[i, col]):
                assert pd.isna(df_out.loc[i, col])


def test_csv_column_order_preservation(desfact_executable, temp_dir):
    """Test that column order is preserved in output CSV."""
    csv_path = Path(temp_dir) / "column_order.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with specific column order
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Property2", "ID", "SMILES", "Property1"])  # Unusual order
        writer.writerow(["b", "1", "CCO", "1.0"])
        writer.writerow(["d", "2", "CC(=O)O", "2.0"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Read original file manually to get exact column order
    with open(csv_path, 'r') as f:
        original_headers = next(csv.reader(f))
    
    # Read output file manually to get column order
    with open(output_path, 'r') as f:
        output_headers = next(csv.reader(f))
    
    # Check that original columns appear in the same order in output
    original_indices = [output_headers.index(col) for col in original_headers]
    assert original_indices == sorted(original_indices)


def test_csv_data_types_preservation(desfact_executable, temp_dir):
    """Test that data types are correctly preserved."""
    csv_path = Path(temp_dir) / "data_types.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with various data types
    data = [
        ["ID", "SMILES", "Integer", "Decimal", "Scientific", "Bool", "Text"],
        ["1", "CCO", "123", "45.67", "1.23e-4", "TRUE", "abc"],
        ["2", "CC(=O)O", "-456", "0.0001", "6.022e23", "FALSE", "def"]
    ]
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Read both CSVs with appropriate type inference
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    # Check data type preservation (as strings to ensure exact format preservation)
    for col in ["Integer", "Decimal", "Scientific", "Bool", "Text"]:
        for i in range(len(df_in)):
            assert str(df_in.loc[i, col]) == str(df_out.loc[i, col])


def test_csv_duplicate_column_handling(desfact_executable, temp_dir):
    """Test handling of duplicate column names in input CSV."""
    csv_path = Path(temp_dir) / "duplicate_columns.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Manually create CSV with duplicate column names
    with open(csv_path, 'w', newline='') as f:
        f.write("ID,SMILES,Property,Property\n")
        f.write("1,CCO,Value1,Value2\n")
        f.write("2,CC(=O)O,Value3,Value4\n")
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # The program should either handle duplicate columns or exit with an error
    if result.returncode == 0:
        # If it continues, check how duplicate columns are handled
        df_out = pd.read_csv(output_path)
        
        # Count columns named "Property" or "Property.X"
        property_cols = [col for col in df_out.columns if col == "Property" or col.startswith("Property.")]
        assert len(property_cols) >= 2  # Should have at least 2 property columns
    else:
        # If it exits with error, check for appropriate message
        assert "duplicate" in result.stderr.lower() or "duplicate" in result.stdout.lower()


def test_csv_recovery_from_partial_processing(desfact_executable, temp_dir):
    """Test recovery behavior when some rows cause processing errors."""
    csv_path = Path(temp_dir) / "partial_errors.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with problematic rows mixed with good rows
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Name"])
        writer.writerow(["1", "CCO", "Good 1"])            # Good
        writer.writerow(["2", "Invalid", "Bad SMILES"])    # Bad
        writer.writerow(["3", "c1ccccc1", "Good 2"])       # Good
        writer.writerow(["4", "", "Empty SMILES"])         # Empty SMILES
        writer.writerow(["5", "CC(=O)O", "Good 3"])        # Good
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # If it processes some rows despite errors:
    if result.returncode == 0 and os.path.exists(output_path):
        df_in = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_path)
        
        # Check if all rows are present
        assert len(df_in) == len(df_out)
        
        # Check if good rows have descriptor values
        good_indices = [0, 2, 4]  # Indices of rows with valid SMILES
        descriptor_cols = [col for col in df_out.columns if col not in df_in.columns]
        
        if len(descriptor_cols) > 0:
            # Check that at least one good row has non-empty descriptor values
            has_descriptor_values = False
            for i in good_indices:
                for col in descriptor_cols:
                    if not pd.isna(df_out.loc[i, col]):
                        has_descriptor_values = True
                        break
                if has_descriptor_values:
                    break
            
            assert has_descriptor_values, "No descriptor values found for valid SMILES"


def test_csv_integrity_with_large_values(desfact_executable, temp_dir):
    """Test CSV integrity with very large text values."""
    csv_path = Path(temp_dir) / "large_values.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with large text values
    large_text = "x" * 10000  # 10,000 character string
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "LargeText"])
        writer.writerow(["1", "CCO", large_text])
        writer.writerow(["2", "CC(=O)O", "Normal text"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check that large text is preserved
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    assert len(df_in.loc[0, "LargeText"]) == len(df_out.loc[0, "LargeText"])
    assert df_in.loc[0, "LargeText"] == df_out.loc[0, "LargeText"]


def test_csv_header_case_sensitivity(desfact_executable, temp_dir):
    """Test handling of case sensitivity in column headers."""
    csv_path = Path(temp_dir) / "case_sensitivity.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with headers that differ only in case
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "property", "Property", "PROPERTY"])
        writer.writerow(["1", "CCO", "lowercase", "Titlecase", "UPPERCASE"])
        writer.writerow(["2", "CC(=O)O", "a", "b", "c"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check if case-differing headers are preserved
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    case_variants = ["property", "Property", "PROPERTY"]
    for header in case_variants:
        assert header in df_out.columns, f"Column {header} not preserved in output"


def test_csv_whitespace_handling(desfact_executable, temp_dir):
    """Test handling of whitespace in CSV fields."""
    csv_path = Path(temp_dir) / "whitespace.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with whitespace in various fields
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Leading", "Trailing", "Both", "Multiple"])
        writer.writerow(["1", "CCO", "   Leading", "Trailing   ", "   Both   ", "Multiple    Spaces"])
        writer.writerow(["2", "CC(=O)O", "\tTab", "Tab\t", "\tTab\t", "Multiple\tTabs"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check if whitespace is preserved or trimmed
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    for i in range(len(df_in)):
        whitespace_cols = ["Leading", "Trailing", "Both", "Multiple"]
        for col in whitespace_cols:
            # Just check that the actual text content is preserved, even if whitespace might be trimmed
            assert df_in.loc[i, col].strip() in df_out.loc[i, col]


def test_csv_comments_in_data(desfact_executable, temp_dir):
    """Test handling of comment-like characters in data fields."""
    csv_path = Path(temp_dir) / "comments.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with comment characters in data
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Comment"])
        writer.writerow(["1", "CCO", "# This looks like a comment"])
        writer.writerow(["2", "CC(=O)O", "// C++ style comment"])
        writer.writerow(["3", "c1ccccc1", "/* Block comment */"])
        writer.writerow(["4", "C1CCCCC1", "-- SQL comment"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check if comment-like strings are preserved
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    comment_styles = ["# This", "// C++", "/* Block", "-- SQL"]
    for i, comment in enumerate(comment_styles):
        assert comment in df_out.loc[i, "Comment"]


def test_csv_handling_of_null_bytes(desfact_executable, temp_dir):
    """Test handling of null bytes in CSV data."""
    csv_path = Path(temp_dir) / "null_bytes.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with null bytes in some fields
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Data"])
        writer.writerow(["1", "CCO", "Normal text"])
        writer.writerow(["2", "CC(=O)O", "Text with \0 null byte"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Program should either handle null bytes or exit with error
    if result.returncode == 0:
        # If it continues, check how null bytes are handled
        with open(output_path, 'r', encoding='utf-8') as f:
            output_content = f.read()
        
        assert "Normal text" in output_content
        # Don't check specifically for null bytes, just ensure the file was processed
    else:
        # If it exits with error, check for appropriate message
        assert any(term in result.stderr.lower() for term in ["null", "invalid", "error"])


def test_csv_handling_of_control_characters(desfact_executable, temp_dir):
    """Test handling of control characters in CSV data."""
    csv_path = Path(temp_dir) / "control_chars.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with control characters
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(["ID", "SMILES", "Control"])
        writer.writerow(["1", "CCO", "Bell: \a"])
        writer.writerow(["2", "CC(=O)O", "Backspace: \b"])
        writer.writerow(["3", "c1ccccc1", "Tab: \t"])
        writer.writerow(["4", "C1CCCCC1", "Vertical tab: \v"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check that non-problematic control chars like tabs are preserved
    df_in = pd.read_csv(csv_path, quoting=csv.QUOTE_ALL)
    df_out = pd.read_csv(output_path, quoting=csv.QUOTE_ALL)
    
    # At minimum, we should have the same number of rows
    assert len(df_in) == len(df_out)
    
    # Check if the control character field names are preserved
    for i in range(len(df_in)):
        assert df_in.loc[i, "Control"].split(":")[0] == df_out.loc[i, "Control"].split(":")[0]


def test_csv_international_delimiter_handling(desfact_executable, temp_dir):
    """Test handling of non-standard delimiters commonly used internationally."""
    # Test different international delimiters
    delimiters = [";", "\t", "|"]
    
    for delimiter in delimiters:
        csv_path = Path(temp_dir) / f"delimiter_{ord(delimiter)}.csv"
        output_path = Path(temp_dir) / f"output_{ord(delimiter)}.csv"
        
        # Create CSV with specific delimiter
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(["ID", "SMILES", "Name"])
            writer.writerow(["1", "CCO", "Ethanol"])
            writer.writerow(["2", "CC(=O)O", "Acetic acid"])
        
        cmd = [
            desfact_executable,
            "-i", str(csv_path),
            "-o", str(output_path),
            "-d", "all",
            "--delimiter", delimiter
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        
        # Check content is preserved
        df_in = pd.read_csv(csv_path, delimiter=delimiter)
        df_out = pd.read_csv(output_path, delimiter=delimiter)
        
        assert len(df_in) == len(df_out)
        assert all(col in df_out.columns for col in df_in.columns)


def test_csv_line_ending_compatibility(desfact_executable, temp_dir):
    """Test compatibility with different line endings (CRLF, LF, CR)."""
    line_endings = {
        "crlf": "\r\n",
        "lf": "\n",
        "cr": "\r"
    }
    
    for name, ending in line_endings.items():
        csv_path = Path(temp_dir) / f"{name}_ending.csv"
        output_path = Path(temp_dir) / f"output_{name}.csv"
        
        # Create CSV with specific line ending
        with open(csv_path, 'w', newline='') as f:
            content = "ID,SMILES,Name" + ending
            content += "1,CCO,Ethanol" + ending
            content += "2,CC(=O)O,Acetic acid" + ending
            f.write(content)
        
        cmd = [
            desfact_executable,
            "-i", str(csv_path),
            "-o", str(output_path),
            "-d", "all"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        
        # Verify that file was processed correctly
        df_in = pd.read_csv(csv_path)
        df_out = pd.read_csv(output_path)
        
        assert len(df_in) == len(df_out)
        assert all(col in df_out.columns for col in df_in.columns)


def test_csv_bom_handling(desfact_executable, temp_dir):
    """Test handling of Byte Order Mark (BOM) in CSV files."""
    csv_path = Path(temp_dir) / "bom.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with BOM
    with open(csv_path, 'wb') as f:
        # Write UTF-8 BOM
        f.write(b'\xef\xbb\xbf')
        
        # Write CSV content
        content = "ID,SMILES,Name\n1,CCO,Ethanol\n2,CC(=O)O,Acetic acid\n"
        f.write(content.encode('utf-8'))
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Verify that file was processed correctly despite BOM
    df_out = pd.read_csv(output_path)
    
    assert len(df_out) == 2
    assert "ID" in df_out.columns
    assert "SMILES" in df_out.columns
    assert "Name" in df_out.columns


def test_csv_handling_of_numeric_precision(desfact_executable, temp_dir):
    """Test preservation of numeric precision for floating point values."""
    csv_path = Path(temp_dir) / "precision.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create CSV with high-precision and special floating point values
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Float", "HighPrecision", "Special"])
        writer.writerow(["1", "CCO", "3.14159", "0.12345678901234567890", "1.0e-20"])
        writer.writerow(["2", "CC(=O)O", "-2.71828", "123456789.123456789", "9.9e+37"])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check precision preservation
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    for i in range(len(df_in)):
        for col in ["Float", "HighPrecision", "Special"]:
            # Compare string representations to check exact format preservation
            assert str(df_in.loc[i, col]) == str(df_out.loc[i, col])


def test_csv_handling_of_quoted_empty_fields(desfact_executable, temp_dir):
    """Test handling of explicitly quoted empty fields versus unquoted empty fields."""
    csv_path = Path(temp_dir) / "quoted_empty.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Manually create CSV with mix of quoted and unquoted empty fields
    with open(csv_path, 'w', newline='') as f:
        f.write('ID,SMILES,Empty1,Empty2\n')
        f.write('1,CCO,,\"\"\n')  # Empty1 is unquoted empty, Empty2 is quoted empty
        f.write('2,CC(=O)O,\"\",\n')  # Empty1 is quoted empty, Empty2 is unquoted empty
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check that all empty fields are treated as empty regardless of quoting
    df_out = pd.read_csv(output_path)
    
    for i in range(2):
        assert pd.isna(df_out.loc[i, "Empty1"]) or df_out.loc[i, "Empty1"] == ""
        assert pd.isna(df_out.loc[i, "Empty2"]) or df_out.loc[i, "Empty2"] == ""


def test_csv_large_descriptor_output(desfact_executable, temp_dir):
    """Test CSV integrity with a large number of descriptor columns."""
    csv_path = Path(temp_dir) / "many_compounds.csv"
    output_path = Path(temp_dir) / "output.csv"
    
    # Create a CSV with 25 compounds (should generate many descriptor columns)
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "SMILES", "Name"])
        
        # Common drug/compound SMILES
        compounds = [
            ("CCO", "Ethanol"),
            ("CC(=O)O", "Acetic acid"),
            ("c1ccccc1", "Benzene"),
            ("C1CCCCC1", "Cyclohexane"),
            ("CC(C)C", "Isobutane"),
            ("COC", "Dimethyl ether"),
            ("CCN", "Ethylamine"),
            ("CC#N", "Acetonitrile"),
            ("c1ccccc1O", "Phenol"),
            ("CC(=O)", "Acetaldehyde"),
            ("CCOC(=O)C", "Ethyl acetate"),
            ("CC(=O)N", "Acetamide"),
            ("CCCC", "Butane"),
            ("c1ccc(cc1)c2ccccc2", "Biphenyl"),
            ("c1ccc2c(c1)cccc2", "Naphthalene"),
            ("c1ccc2c(c1)ccc3c2cccc3", "Anthracene"),
            ("Cc1ccccc1", "Toluene"),
            ("ClC(Cl)(Cl)Cl", "Carbon tetrachloride"),
            ("CC(=O)OC1C(C(C(C(O1)OC2C(C(C(C(O2)OC3C(C(C(C(O3)C)(C)O)O)O)(C)O)O)O)(C)O)O)O", "Vancomycin"),
            ("CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC", "Cocaine"),
            ("C1C(C(OC2=CC=CC=C21)CO)O", "Epinephrine"),
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine"),
            ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "Ibuprofen"),
            ("CC12CCC(CC1)C(C)(C)O2", "Eucalyptol"),
            ("OC1=C(C=CC(=C1)O)O", "Pyrogallol")
        ]
        
        for i, (smiles, name) in enumerate(compounds, 1):
            writer.writerow([str(i), smiles, name])
    
    cmd = [
        desfact_executable,
        "-i", str(csv_path),
        "-o", str(output_path),
        "-d", "all"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    
    # Check output has all molecules and many descriptor columns
    df_in = pd.read_csv(csv_path)
    df_out = pd.read_csv(output_path)
    
    assert len(df_in) == len(df_out)
    assert all(col in df_out.columns for col in df_in.columns)
    
    # Should have significantly more columns in output
    descriptor_cols = [col for col in df_out.columns if col not in df_in.columns]
    assert len(descriptor_cols) > 5, "Not enough descriptor columns generated"


def test_csv_nonexistent_descriptor(desfact_executable, valid_smiles_csv, temp_dir):
    """Test behavior when requesting a non-existent descriptor."""
    output_path = Path(temp_dir) / "output.csv"
    
    cmd = [
        desfact_executable,
        "-i", str(valid_smiles_csv),
        "-o", str(output_path),
        "-d", "nonexistent_descriptor"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Should either exit with error or continue with warning
    if result.returncode != 0:
        assert "unknown" in result.stderr.lower() or "invalid" in result.stderr.lower() or "not found" in result.stderr.lower()
    else:
        # If it continues, check warning message and output integrity
        assert "unknown" in result.stdout.lower() or "unknown" in result.stderr.lower()
        
        df_in = pd.read_csv(valid_smiles_csv)
        df_out = pd.read_csv(output_path)
        
        # Output should still preserve original data
        assert len(df_in) == len(df_out)
        assert all(col in df_out.columns for col in df_in.columns)