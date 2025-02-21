"""
Given a sequence of a CAS9 target site, this script will generate the oligos incorporating Hammerhead to be cloned in pNOC-Cas9 vector.
Usage: python generate_CAS9_primers.py <target_sequence> <case_formatting>
Arguments:
    target_sequence: The target sequence of the CAS9 site
    case_formatting: True if you want the oligos to swap between upper and lower case, False otherwise
Optional:
    -h, --help: Show this help message and exit
"""

def help():
    print(__doc__)

from guide_helper import ordering_info
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        help()
        sys.exit(1)
    target = sys.argv[1]
    case_formatting = False if len(sys.argv) < 3 else sys.argv[2].lower() == "true"
    print( ordering_info(target, case_formatting=case_formatting ) )