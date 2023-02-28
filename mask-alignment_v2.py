"""
Mask initial bases from alignment FASTA
"""

# this script was originaly created by nextstrain and distributed at https://github.com/nextstrain/ncov/blob/master/scripts/mask-alignment.py
# this version is modified in that it takes as input a file listing sites to mask and not a list in plain text


import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument("--mask-from-beginning", type = int, required=True, help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type = int, help="number of bases to mask from end")
    parser.add_argument("--mask-sites", type = str,  help="list of sites to mask")
    parser.add_argument("--output", required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    being_length = 0
    if args.mask_from_beginning:
        begin_length = args.mask_from_beginning
    end_length = 0
    if args.mask_from_end:
        end_length = args.mask_from_end

    f = open(args.mask_sites, "r")
    arr_intValues = []
    for myLine in f:
     arr_intValues.append(int(myLine))

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            start = "N" * begin_length
            middle = seq[begin_length:-end_length]
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if args.mask_sites:
                for site in arr_intValues:
                    seq_list[site-1] = "N"
            record.seq = Seq("".join(seq_list))
            Bio.SeqIO.write(record, outfile, 'fasta')
