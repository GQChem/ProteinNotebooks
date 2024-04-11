import sys

def find_cpio_offset(file_content):
    """
    Finds the offset of the CPIO archive in the RPM file by searching for
    the magic numbers '070701' or '070702' (hex representation of CPIO newc format).
    Returns the offset, or -1 if not found.
    """
    # Searching for the magic numbers indicating the start of a CPIO archive
    offset = file_content.find(b'070701')
    if offset == -1:
        offset = file_content.find(b'070702')
    return offset

def extract_cpio_from_rpm(file_path):
    """
    Extracts the CPIO archive from an RPM file, handling '-' as stdin.
    Outputs the CPIO archive to stdout.
    """
    # Determine whether to read from stdin or a file
    if file_path == '-':
        # Read all bytes from stdin
        rpm_content = sys.stdin.buffer.read()
    else:
        # Open the file and read its content
        with open(file_path, 'rb') as rpm_file:
            rpm_content = rpm_file.read()

    # Find the offset of the CPIO archive
    offset = find_cpio_offset(rpm_content)
    if offset == -1:
        raise ValueError("CPIO start not found in RPM file.")

    # Output the CPIO archive to stdout
    sys.stdout.buffer.write(rpm_content[offset:])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rpm2cpio.py <rpm_file or ->", file=sys.stderr)
        sys.exit(1)

    try:
        extract_cpio_from_rpm(sys.argv[1])
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
