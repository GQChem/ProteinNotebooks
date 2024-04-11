import sys

def find_cpio_offset(rpm_path):
    with open(rpm_path, 'rb') as rpm_file:
        rpm_content = rpm_file.read()
    # The cpio archive starts with '070701' or '070702' (for new cpio format)
    offset = rpm_content.find(b'070701') 
    if offset == -1:
        offset = rpm_content.find(b'070702')
    return offset

def extract_cpio(rpm_path, offset):
    with open(rpm_path, 'rb') as rpm_file:
        rpm_file.seek(offset)
        while True:
            chunk = rpm_file.read(4096)
            if not chunk:
                break
            sys.stdout.buffer.write(chunk)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rpm2cpio.py <rpm_file>")
        sys.exit(1)
    rpm_path = sys.argv[1]
    offset = find_cpio_offset(rpm_path)
    if offset != -1:
        extract_cpio(rpm_path, offset)
    else:
        print("Error: CPIO start not found.")
        sys.exit(2)