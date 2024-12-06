import re


def calculate_identity(cigar_string):
    """
    Calculate the identity of a sequence based on its CIGAR string.

    Args:
        cigar_string (str): CIGAR string (e.g., "17S3=1X21=1D8=1D23=1I").

    Returns:
        float: Identity value as a ratio.
    """
    # Parse the CIGAR string using regex
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)

    # Initialize counts
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0

    # Process each operation
    for length, op in operations:
        length = int(length)
        if op == "=":
            matches += length
        elif op == "X":
            mismatches += length
        elif op == "I":
            insertions += length
        elif op == "D":
            deletions += length

    # Compute identity
    total_bases = matches + mismatches + insertions + deletions
    if total_bases == 0:
        return 0  # Avoid division by zero

    identity = matches / total_bases
    return identity


def calculate_reference_length(cigar_string):
    """
    Calculate the reference length of a sequence alignment based on its CIGAR string.

    Args:
        cigar_string (str): CIGAR string (e.g., "17S3=1X21=1D8=1D23=1I").

    Returns:
        int: Reference length.
    """
    # Parse the CIGAR string using regex
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)

    # Initialize reference length
    reference_length = 0

    # Process each operation
    for length, op in operations:
        length = int(length)
        if op in ("=", "X", "D"):  # Contributes to reference length
            reference_length += length

    return reference_length


if __name__ == "__main__":

    # 示例调用
    cigar = "17=1I3=2I18=1I4=1D20=4I5=3I15=1I6=1I4=1D22=1I5=1D16=1D7=1I8=2X5=1I48=2I9=1I2=1D20=1I12=1D10=1X30=1D4=1D9=1I42=1I4=1D27=2D3=1D21=1D3=1D43=1D7=2D11=1D14=2I8=1D8=2X15=2D5=2I3=1D24=1D18=1D3=1D7=1D21=1X10=1D5=1I36=1I28=1X1=2D2=1D1=1D5=1I14=1I21=1I5=1D5=1I10=1D3=2D24=1I9=1D2=2D12=1D15=1D16=1D3=1D3=1D4=1D2=1D22=1I13=1I30=1D11=1I21=1D15=1I6=1I9=1D32=2D4=1X1=1X5=1I26=3I11=1I5=1D27=2D12=2I6=1D5=1I11=1D2=1D9=1D7=1I3=1I18=1D16=1X1=1X7=1I39=1X4=1D5=2I18=1D6=1D19=1I16=1I3=1X5=2I5=2X32=1I8=1D10=1D6=1D16=1I35=1D4=1I4=1I17=1I34=1D8=1D21=2D5=1D8=1D21=1I15=2I14=1D14=1X2=1D6=2I5=1I22=1D4=1D3=1I23=2D8=2I11=1X7=1D15=2I2=1D10=1D9=2I13=1I13=4I2=10I2=5I13=1I11=8S"
    identity = calculate_identity(cigar)
    print(f"CIGAR: {cigar}")
    print(f"Identity: {identity:.2%}")
    length = calculate_reference_length(cigar)
    print(f"ref_aligned_length:{length}")

    cigar = "17=1I3=2I18=1I4=1D20=4I5=3I15=1I6=1I4=1D22=1I5=1D16=1D7=1I8=2X5=1I48=2I9=1I2=1D20=1I12=1D10=1X30=1D4=1D9=1I42=1I4=1D27=2D3=1D21=1D3=1D43=1D7=2D11=1D14=2I8=1D8=2X15=2D5=2I3=1D24=1D18=1D3=1D7=1D21=1X10=1D5=1I36=1I28=1X1=2D2=1D1=1D5=1I14=1I21=1I5=1D5=1I10=1D3=2D24=1I9=1D2=2D12=1D15=1D16=1D3=1D3=1D4=1D2=1D22=1I13=1I30=1D11=1I21=1D15=1I6=1I9=1D32=2D4=1X1=1X5=1I26=3I11=1I5=1D27=2D12=2I6=1D5=1I11=1D2=1D9=1D7=1I3=1I18=1D16=1X1=1X7=1I39=1X4=1D5=2I18=1D6=1D19=1I16=1I3=1X5=2I5=2X32=1I8=1D10=1D6=1D16=1I35=1D4=1I4=1I17=1I34=1D8=1D21=2D5=1D8=1D21=1I15=2I14=1D14=1X2=1D6=2I5=1I22=1D4=1D3=1I23=2D8=2I11=1X7=1D15=2I2=1D10=1D9=2I13=1I13=56S"
    identity = calculate_identity(cigar)
    print(f"CIGAR: {cigar}")
    print(f"Identity: {identity:.2%}")
    length = calculate_reference_length(cigar)
    print(f"ref_aligned_length:{length}")
