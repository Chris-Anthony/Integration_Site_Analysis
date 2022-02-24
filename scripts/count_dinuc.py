#!/usr/bin/env python

print("Initializing...")
print()

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes a fasta file and counts the begining dinucleotides.""")
    parser.add_argument('-f','--file',  action='store', nargs='?', type=str, 
                    required=True, help='fasta file of mapped reads', dest="fh")

    return parser.parse_args()

def main(fh:str):
    seq = "NN"

    dinuc_dict = dict()

    with open(fh, "rt") as fh:
        for line in fh:
            if line[0] != ">":
                seq += line.strip()
            else:
                if seq[0:2] in dinuc_dict:
                    dinuc_dict[seq[0:2]] += 1
                else: 
                    dinuc_dict[seq[0:2]] = 1

                seq = str()
    
    print("test")

    for key in dinuc_dict:
        print(key, "\t", dinuc_dict[key])

if __name__ =="__main__":
    args = get_args()
    fh = args.fh

    main(fh)