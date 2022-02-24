#!/usr/bin/env python

print("Initializing...")

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes a bedgraph file and counts number of integration sites.""")
    parser.add_argument('-x', '--input',  action='store', nargs='?', type=str, 
                    required=True, help='bedgraph file that ends in .bedgraph', dest="bg")
    parser.add_argument('-c', '--coverage-cutoff',  action='store', nargs='?', type=int, 
                    required=True, help='ignore locations with less than this coverage cutoff', dest="c")
    parser.add_argument('-l', '--length-cutoff',  action='store', nargs='?', type=int, 
                    required=True, help='ignore locations with lengths less than this cutoff', dest="l")
    parser.add_argument('-o', '--outdir',  action='store', nargs='?', type=str, 
                    required=False, help='output directory', dest="outdir")
    parser.add_argument('-s', '--strand',  action='store', nargs='?', type=str, default="both", 
                    required=False, help='strandedness', dest="strand")

    return parser.parse_args()

def add_site(sites:dict, chr:str, start:str, end:str, sum_cov:int, length:int, dir:str)->dict:
    if length == 0:
        length = 1

    if chr in sites:
        sites[chr].append((int(start), int(end), sum_cov/length, length, dir))
    else:
        sites[chr]=list()
        sites[chr].append((int(start), int(end), sum_cov/length, length, dir))

    return sites


def main(bg:str, c:int, l:int, outdir:str, strand:str):
    import re

    # obtain file name and open output files
    name = re.match("[.\/\w]+/(.+)\.bedgraph", bg).group(1)
    outfile_IS = open(outdir + "/" + name + "_IS" + ".txt", "w")
    outfile_coor_fwd = open(outdir + "/" + name + "_fwd.bed", "w")
    outfile_coor_rev = open(outdir + "/" + name + "_rev.bed", "w")

    sites = dict()

    chr, dir = str(), str()
    start, end = 0, 0
    cov_start, cov_end, sum_cov = 1, 1, 0
    counter_start, counter_max = 0, 0

    # internal breakpoints if integration sites overlap
    min_cov = [1000000000, 0, 0]  # value, coordinate
    max_cov = [0, 0, 0]  # value, coordinate

    with open(bg, "rt") as fh:
        for line in fh:
            line_list = line.split("\t")

            if end != line_list[1]: # if the end of the previous line is NOT equal to current line
                # print(start)

                if strand == "both":
                    if (0.9*max_cov[0]) > cov_start and (0.9*max_cov[0]) > cov_end: # shape "^"
                        # print("I'm doing this - ^")
                        # do the first half
                        dir = "-"
                        sites = add_site(sites, chr, start, max_cov[1], max_cov[2], int(max_cov[1])-int(start), dir)

                        # do the second half
                        dir = "+"
                        sites = add_site(sites, chr, max_cov[1], end, sum_cov-max_cov[2], int(end)-int(max_cov[1]), dir)

                    elif (1.1*min_cov[0]) < cov_start and (1.1*min_cov[0]) < cov_end: # shape "V"
                        # print("I'm doing this - V")
                        # do the 1st half
                        dir = "+"
                        sites = add_site(sites, chr, start, min_cov[1], min_cov[2], int(min_cov[1])-int(start), dir)

                        # do the second half
                        dir = "-"
                        sites = add_site(sites, chr, min_cov[1], end, sum_cov-min_cov[2], int(end)-int(min_cov[1]), dir)

                    else:
                        if min_cov[1] < max_cov[1]:
                            # print("I'm doing this - /")
                            dir = "-"
                        else: # cov_start > cov_end:
                            # print("I'm doing this - \\")
                            dir = "+"
                        
                        sites = add_site(sites, chr, start, end, sum_cov, int(end)-int(start), dir)
                elif strand == "+":               
                    dir = "+"
                    sites = add_site(sites, chr, start, end, sum_cov, int(end)-int(start), dir)

                elif strand == "-":
                    dir = "-"
                    sites = add_site(sites, chr, start, end, sum_cov, int(end)-int(start), dir)

                chr = line_list[0]
                start = line_list[1]

                # reset
                cov_start = int(float(line_list[3].rstrip()))
                cov_end = 1
                sum_cov = 0

                min_cov = [1000000000, 0, 0]
                max_cov = [0, 0, 0]

                counter_max = 0
                counter_start = 0

            if counter_start <= 3: # accounts for misaligned reads at the 5' end
                if cov_start/int(float(line_list[3].rstrip())) < 0.01:
                    cov_start = int(float(line_list[3].rstrip()))
                    start = line_list[1]

                    min_cov[0] = int(float(line_list[3].rstrip()))
                    min_cov[1] = end

                counter_start += 1

            if int(float(line_list[3].rstrip()))/cov_end > 0.01: # accounts for misaligned reads at 3' end
                cov_end = int(float(line_list[3].rstrip()))
                end = line_list[2]

                sum_cov += (cov_end)*(int(end)-int(line_list[1]))

            if (1.1*min_cov[0]) > cov_end:
                min_cov[0] = cov_end
                min_cov[1] = end
                min_cov[2] += (cov_end)*(int(end)-int(line_list[1]))

            if (0.9*max_cov[0]) < cov_end:
                if (0.9*max_cov[0]) < cov_end and cov_end < (1.1*max_cov[0]): # coverage has stabilized
                    counter_max += 1
                
                if counter_max < 3:
                    max_cov[0] = cov_end
                    max_cov[1] = end
                    max_cov[2] += (cov_end)*(int(end)-int(line_list[1]))

                else: # counter > 3
                    if (0.9*max_cov[0]) < cov_end < (1.1*max_cov[0]):  # handles both scenarios where right or left edge is bigger
                        max_cov[0] = cov_end
                        max_cov[1] = end
                        max_cov[2] += (cov_end)*(int(end)-int(line_list[1]))

    counts = dict()

    for key in sites:
        if key != "":
            for IS in sites[key]:
                if IS[2] >= c:
                    if IS[3] >= l:
                        if str(IS[4]) == "+":
                            outfile_coor_fwd.write(key + "\t" + str(IS[0]) + "\t" + str(IS[1]) + "\t" + "\t" + "\t" + "+" + "\n")
                        else: # if direction is "-"
                            outfile_coor_rev.write(key + "\t" + str(IS[0]) + "\t" + str(IS[1]) + "\t" + "\t" + "\t" + "-" + "\n")
                        
                        if key in counts:
                            counts[key] += 1
                        else:
                            counts[key] = 1

    outfile_IS.write("Number of Integration Sites Per Chromosome\n")
    for key in counts:
        if key != "":
            outfile_IS.write(key + ":\t" + str(counts[key]) + "\n")

    outfile_IS.write("\nGenomic Locations of Integration Sites\n")
    outfile_IS.write("chr\tstart\tend\tdirection\tavg.cov\tlength\n")
    for key in sites:
        if key != "":
            for IS in sites[key]:
                if IS[2] >= c:
                    if IS[3] >= l:
                        outfile_IS.write(key + "\t" + str(IS[0]+1) + "\t" + str(IS[1]+1) + "\t" + str(IS[4]) + "\t" + str("{:.2f}".format(IS[2])) + "\t" + str(IS[3]) + "\n" )


    outfile_IS.close()
    outfile_coor_fwd.close()
    outfile_coor_rev.close()

    print("Complete.")


if __name__ =="__main__":
    args = get_args()
    bg, c, l, outdir, strand = args.bg, args.c, args.l, args.outdir, args.strand

    main(bg, c, l, outdir, strand)
