#!/usr/bin/env python
print("Initializing...")

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes a bedgraph file and counts number of integration sites.""")
    parser.add_argument('-i', '--input',  action='store', nargs='?', type=str, 
                    required=True, help='sam file', dest="sam")
    parser.add_argument('-c', '--coverage-cutoff',  action='store', nargs='?', type=int, default=5,
                    required=False, help='ignore locations with less than this coverage cutoff', dest="c")
    parser.add_argument('-l', '--local',  action='store', nargs='?', type=int, default=10,
                    required=False, help='how near the next read has to be in order to be part of group', dest="l")
#    parser.add_argument('-l', '--length-cutoff',  action='store', nargs='?', type=int, 
#                    required=False, help='ignore locations with lengths less than this cutoff', dest="l")
    parser.add_argument('-o', '--outdir',  action='store', nargs='?', type=str, default="./",
                    required=False, help='output directory', dest="outdir")

    return parser.parse_args()

def main(sam:str, c:int, l:int, outdir:str):
    import re

    # obtain file name and open output files
    name = re.match("[.\/\w]+/(.+)\.sam", sam).group(1)
    outfile_coor_fwd = open(outdir + "/" + name + "_fwd.bed", "w")
    outfile_coor_rev = open(outdir + "/" + name + "_rev.bed", "w")

    # initalize data structures
    pos = 0

    dict_rev_sites = dict() # max_position are keys, values are dictionary, # possible issues keys will be constantly updating
    dict_rev_sites[pos] = dict()
    dict_rev_sites[pos][pos] = dict()
    dict_rev_sites[pos][pos][0] = 0

    fwd_sites = dict()
    fwd_sites[pos] = dict()
    fwd_sites[pos][0] = 0

    # initialize temp variables
    rev_chr = str()

    fwd_chr = str()
    fwd_pos = int()

    print("Reading...")

    with open(sam, "rt") as fh:
        for line in fh:
            if line[0] != "@":
                line_list = line.split("\t")
                flag = int(line_list[1])
                new_chr = str(line_list[2])
                new_pos = int(line_list[3])
                length = len(line_list[9])

                if flag & 16 == 16: # for reverse strands
                    new_pos = new_pos + length

                    if new_chr == rev_chr:
                        in_list = False

                        # check if the new position is within l bases of the max positions within dictionary
                        for temp_pos in dict_rev_sites:
                            if temp_pos - l <= new_pos <= temp_pos + l:
                                in_list = True
                                max_pos = temp_pos
                        
                        if in_list:
                            # check if new position exists within max_pos dictionary
                            if new_pos in dict_rev_sites[max_pos]:
                                pass
                            else:
                                dict_rev_sites[max_pos][new_pos] = dict()

                            # add count to length within pos dictionary
                            if length in dict_rev_sites[max_pos][new_pos]:
                                dict_rev_sites[max_pos][new_pos][length] += 1
                            else:
                                dict_rev_sites[max_pos][new_pos][length] = 1

                            # update the max_pos 
                            max_count = 0
                            temp_pos = 0

                            for site in dict_rev_sites[max_pos]:
                                for length in dict_rev_sites[max_pos][site]:
                                    if dict_rev_sites[max_pos][site][length] > max_count:
                                        max_count = dict_rev_sites[max_pos][site][length]
                                        temp_pos = site
                            
                            dict_rev_sites[temp_pos] = dict_rev_sites.pop(max_pos)

                        else:
                            if len(dict_rev_sites) == 5:
                                # export data
                                for max_pos in dict_rev_sites:
                                    i = 0

                                    for site in dict_rev_sites[max_pos]:
                                        j = 0

                                        for k in list(dict_rev_sites[max_pos][site].values()):
                                            j += k

                                        if j > i:
                                            pos = site
                                            i = j
                                        elif j == i:
                                            pos = site
                                    
                                    key_list = list(dict_rev_sites[max_pos][pos].keys())
                                    val_list = list(dict_rev_sites[max_pos][pos].values())

                                    if max(val_list) >= c:
                                        seqlen = key_list[val_list.index(max(val_list))]
                                        start = pos - seqlen
                                        
                                        outfile_coor_rev.write(rev_chr + "\t" + str(start) + "\t" + str(pos - 1) + "\t" + "\t" + "\t" + "-" + "\n")

                                #reset
                                dict_rev_sites = dict()
                                dict_rev_sites[new_pos] = dict()
                                dict_rev_sites[new_pos][new_pos] = dict()
                                dict_rev_sites[new_pos][new_pos][length] = 1
                            else:
                                # add new_position to the dictionary
                                dict_rev_sites[new_pos] = dict()
                                dict_rev_sites[new_pos][new_pos] = dict()
                                dict_rev_sites[new_pos][new_pos][length] = 1
                    else:
                        #export
                        for max_pos in dict_rev_sites:
                            i = 0

                            for site in dict_rev_sites[max_pos]:
                                j = 0

                                for k in list(dict_rev_sites[max_pos][site].values()):
                                    j += k

                                if j > i:
                                    pos = site
                                    i = j
                                elif j == i:
                                    pos = site
                            
                            key_list = list(dict_rev_sites[max_pos][pos].keys())
                            val_list = list(dict_rev_sites[max_pos][pos].values())

                            if max(val_list) >= c:
                                seqlen = key_list[val_list.index(max(val_list))]
                                start = pos - seqlen
                                
                                outfile_coor_rev.write(rev_chr + "\t" + str(start) + "\t" + str(pos - 1) + "\t" + "\t" + "\t" + "-" + "\n")

                        #reset
                        dict_rev_sites = dict()
                        dict_rev_sites[new_pos] = dict()
                        dict_rev_sites[new_pos][new_pos] = dict()
                        dict_rev_sites[new_pos][new_pos][length] = 1

                        rev_chr = new_chr

                else: # for forward strands
                    if new_chr != fwd_chr or new_pos > fwd_pos + l:
                        # multiple loci in the group, choose the one with the highest read counts
                        i = 0

                        for site in fwd_sites:
                            j = 0

                            for k in list(fwd_sites[site].values()):
                                j += k

                            if j > i:
                                pos = site
                                i = j
                            elif j == i:
                                pos = site
                        
                        key_list = list(fwd_sites[pos].keys())
                        val_list = list(fwd_sites[pos].values())

                        if max(val_list) >= c:
                            seqlen = key_list[val_list.index(max(val_list))]
                            end = pos + seqlen
                            
                            outfile_coor_fwd.write(fwd_chr + "\t" + str(pos) + "\t" + str(end - 1) + "\t" + "\t" + "\t" + "+" + "\n")

                        #reset
                        fwd_sites = dict()
                        fwd_sites[new_pos] = dict()
                        fwd_sites[new_pos][length] = 1

                    else:
                        if new_pos in fwd_sites:
                            pass
                        else:
                            fwd_sites[new_pos] = dict()

                        if length in fwd_sites[new_pos]:
                            fwd_sites[new_pos][length] += 1
                        else:
                            fwd_sites[new_pos][length] = 1

                    fwd_chr = new_chr
                    fwd_pos = new_pos

    # last reverse site
    for max_pos in dict_rev_sites:
        i = 0

        for site in dict_rev_sites[max_pos]:
            j = 0

            for k in list(dict_rev_sites[max_pos][site].values()):
                j += k

            if j > i:
                pos = site
                i = j
            elif j == i:
                pos = site
        
        key_list = list(dict_rev_sites[max_pos][pos].keys())
        val_list = list(dict_rev_sites[max_pos][pos].values())

        if max(val_list) >= c:
            seqlen = key_list[val_list.index(max(val_list))]
            start = pos - seqlen
            
            outfile_coor_rev.write(rev_chr + "\t" + str(start) + "\t" + str(pos - 1) + "\t" + "\t" + "\t" + "-" + "\n")

    # last foward site
    i = 0

    for site in fwd_sites:
        j = 0

        for k in list(fwd_sites[site].values()):
            j += k

        if j > i:
            pos = site
            i = j
        elif j == i:
            pos = site
    
    key_list = list(fwd_sites[pos].keys())
    val_list = list(fwd_sites[pos].values())

    if max(val_list) >= c:
        seqlen = key_list[val_list.index(max(val_list))]
        end = pos + seqlen
        
        outfile_coor_fwd.write(fwd_chr + "\t" + str(pos) + "\t" + str(end - 1) + "\t" + "\t" + "\t" + "+" + "\n")

    # close files
    outfile_coor_fwd.close()
    outfile_coor_rev.close()

    print("Complete.")

if __name__ =="__main__":
    args = get_args()
    sam, c, l, outdir = args.sam, args.c, args.l, args.outdir

    main(sam, c, l, outdir)
