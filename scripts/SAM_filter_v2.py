#!/usr/bin/env python

print("Initializing...")
print()

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes a SAM file and filters passing reads based on user inputs.""")
    parser.add_argument('-x','--input',  action='store', nargs='?', type=str, 
                    required=True, help='sam file that ends in .sam', dest="sam")
    parser.add_argument('-m','--mapq',  action='store', nargs='?', type=int, 
                    required=True, help='mapq threshold that reads need to be above to pass filter', dest="mapq")
    parser.add_argument('-d','--dinucleotide',  action='store', nargs='?', type=str, 
                    required=True, help='dinucleotide that reads need to begin with to pass filter', dest="dinuc")

    return parser.parse_args()

def main(sam:str, mapq:int, dinuc:str):
    '''Takes in sam file, filters out reads that are below mapq threshold or start with wrong dinucleotide'''
    import re

    low_mapping = 0
    wrong_dinucleotide = 0
    pass_all = 0

    # obtain file name and open output files
    name = re.match("(.+)\.sam", sam).group(1)
    outfile_pass = open(name + "_pass" + ".sam", "w")
    outfile_fail = open(name + "_fail" + ".sam", "w")

    with open(sam, "rt") as fh:
        for line in fh:
            if line[0] == "@":
                outfile_pass.write(line)
                outfile_fail.write(line)
            else:
                line_list = line.split("\t")

                mapping_score = int(line_list[4]) # pull mapq score
                seq = str(line_list[9]) # pull sequence score

                if mapping_score < mapq:
                    low_mapping += 1
                    outfile_fail.write(line)

                else:
                    if dinuc == "no":
                        outfile_pass.write(line)
                        pass_all += 1
                    else:
                        if seq[0:2]!=dinuc:
                            outfile_fail.write(line)
                            wrong_dinucleotide += 1
                        else:
                            outfile_pass.write(line)
                            pass_all += 1

    # Close all files
    outfile_pass.close()
    outfile_fail.close()

    print("---SUMMARY---")
    print("low mapping:\t", low_mapping)
    print("wrong dinucleotide:\t", wrong_dinucleotide)
    print("no. passing reads:\t", pass_all)

if __name__ =="__main__":
    args = get_args()
    sam = args.sam
    mapq = args.mapq
    dinuc = args.dinuc

    main(sam, mapq, dinuc)