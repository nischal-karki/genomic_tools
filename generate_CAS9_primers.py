"""
Given a sequence of a CAS9 target site, this script will generate the oligos incorporating Hammerhead to be cloned in pNOC-Cas9 vector.
Usage: python generate_CAS9_primers.py <target_sequence> <case_formatting>
Arguments:
    target_sequence: The target sequence of the CAS9 site
    case_formatting: True if you want the oligos to swap between upper and lower case, False otherwise
Optional:
    -h, --help: Show this help message and exit
    -c [cas-type], --cas=[cas-type]: Specify the type of CAS9 to use (default: CAS9)
    Supported CAS types:
        SpCas9, SpCas9_VRER, SpCas9_VQR, xCas9, SpCas9_NG, SaCas9, AsCpf1, AsCpf1_RR, AsCpf1_RVR, LbCpf1, LbCpf1_RR, FnCas12a, CjCas9, NmCas9, StCas9, TdCas9
"""

def help():
    print(__doc__)

from guide_helper import oligos_to_order
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        help()
        sys.exit(1)
    target = sys.argv[1]
    case_formatting = False if len(sys.argv) < 3 else sys.argv[2].lower() == "true"
    five, three = oligos_to_order(target, case_formatting)
    print(f"fwd: 5'{five}3'")
    print(f"rev: 5'{three}3'")