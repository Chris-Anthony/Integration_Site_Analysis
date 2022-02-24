#!/usr/bin/env python

print("Starting...")

def get_args():
    '''
    Argparse code to allow for in-line command interface. This script requires pysam.
    '''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes a fasta file and filters seqeunces for start sites.""")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i','--input',  action='store', nargs='?', type=str, 
                    required=True, help='bed file that ends in .bed', dest="inp")
    required.add_argument('-x','--seq',  action='store', nargs='?', type=str, 
                    required=True, help='true start site sequence', dest="seq")
    required.add_argument('-s','--size',  action='store', nargs='?', type=str, 
                    required=True, help='ChromSize file that contains the size of each chromosome', dest="size")
    required.add_argument('-b','--build',  action='store', nargs='?', type=str, 
                    required=True, help='fasta file of the genome build', dest="build")

    optional.add_argument('-o', '--outdir',  action='store', nargs='?', type=str, default="./",
                    required=False, help='output directory', dest="outdir")

    # adds or removes bases up to the threshold from the start site
    optional.add_argument('-t','--threshold',  action='store', nargs='?', type=int, default=5,
                    required=False, help='''Maximum number of bases to check upstream or downstream of the beginning 
                    of the read to search for the start site. The default is 5 bases.''', dest="thres")
    optional.add_argument('-l','--insert-length',  action='store', nargs='?', type=int, default=0,
                    required=False, help='''Theoretical insert length. Modifications that yield a read closer to insert length is prefered''', 
                    dest="length")
    
    # more specific 5' addition or 5' removal from the start site
    optional.add_argument('-u','--upstream',  action='store', nargs='?', type=int, default=0,
                    required=False, help='''Maximum number of bases to check upstream of the beginning 
                    of the read to search for the start site''', dest="up")
    optional.add_argument('-d','--downstream',  action='store', nargs='?', type=int, default=0,
                    required=False, help='''Maximum number of bases to check downstream of the beginning 
                    of the read to search for the start site''', dest="down")

    return parser.parse_args()

def revComp(seq:str)->str:
    '''
    Take in a sequence string. Returns its reverse complement.
    '''
    
    # Reverse the sequence
    rev = seq[::-1]
    
    # Watson-Crick base pair as a dictionary
    convert_dict = {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C",
        "N":"N"
        }

    # Watson-Crick base pair conversion 
    for i in range(len(rev)):
        if i == 0:
            revcompseq = convert_dict[rev[i]]
        else:
            revcompseq = f'{revcompseq}{convert_dict[rev[i]]}'

    return revcompseq

def checkUp(read:str, seq:str, count:int)->int:
    '''
    Takes in a read. Aligns the sequence to the start of the read. Reports any matches to the 3' of the sequence.
    '''
    if read[0:(len(seq)-count)] == seq[count:]:
        pass
    else:
        count += 1
        count = checkUp(read, seq, count)

    return count

def main(inp:str, seq:str, size:str, build:str, outdir:str, thres:int, length:int, up:int, down:int):
    '''
    Takes in fasta file, filters out reads do not start with a specified sequence. 
    Attempts to add bases upstream of read to recover the specified sequence
    '''
    import re
    import pysam

    # open fasta file to obtain sequences 
    openFasta = pysam.Fastafile(build)

    # Get the sizes of the chromosomes
    chromSize = dict()

    with open(size, "rt") as fh:
        for line in fh:
            line = line.split()
            chromSize[line[0].upper().strip("CHR")] = int(line[1])

    # Open output file
    name = re.match("[.\/\w]+/(.+)\.bed", inp).group(1)
    outfile = open(outdir + "/" + name + "_recovered" + ".bed", "w")
    notfile = open(outdir + "/" + name + "_unrecovered" + ".bed", "w")

    # set the maximum number of bases added or removed to the beginning of the read to find the start site
    if up or down:
        pass
    else:
        up = thres
        down = thres

    print("Reading...")

    with open(inp, "rt") as fh:
        for line in fh:
            # get coordinates from line
            getCoors = line.split('\t')
            chr = str(getCoors[0])
            start = int(getCoors[1])
            end = int(getCoors[2])
            strand = str(getCoors[5]).strip()

            # obtain the read
            read = openFasta.fetch(chr, start, end+1) # openFasta is right-end exclusive

            if strand == "-":
                read = revComp(read)

            # initialize
            temp_start = start
            temp_end = end
            temp_read = read 

            # search upstream
            match_up = False
            count_up = 0

            while not match_up:
                # try to align the sequence of interest to the read
                count = checkUp(read, seq, 0)

                if count == 0: # if there are no changes to counts, exit loop
                    match_up = True
                else:
                    if count + count_up > up: # check whether adding count it over the threshold
                        count_up = -999
                        match_up = True
                    else:
                        count_up += count

                # make sure update to coordinates are within the size of the chromosomes
                if strand == "-":
                    if end + count > chromSize[chr]:
                        count_up = -999
                        match_up = True
                    else:
                        end += count
                else:
                    if start - count < 0:
                        count_up = -999
                        match_up = True
                    else:
                        start -= count

                # obtain sequence for new coordinates
                read = openFasta.fetch(chr, start, end+1)

                if strand == "-":
                    read = revComp(read)

            # search internally/ downstream of read start 
            #print(temp_read)
            count_down = re.search(seq, temp_read)

            if count_down == None:
                count_down = 999
            else:
                count_down = count_down.start()

                if count_down > down:
                    count_down = 999
                else:
                    if strand == "-":
                        temp_end -= count_down
                    else:
                        temp_start += count_down

            # write coordinates to output bed file
            if count_up == -999 and count_down == 999:
                notfile.write(f'{chr}\t{temp_start}\t{temp_end}\t{count_up}\t\t{strand}\n')
            elif count_up == -999 and count_down != 999:
                outfile.write(f'{chr}\t{temp_start}\t{temp_end}\tremoved_{count_down}_bases\t\t{strand}\n')
            elif count_up != -999 and count_down == 999:
                outfile.write(f'{chr}\t{start}\t{end}\tadded_{count_up}_bases\t\t{strand}\n')
            else:
                if length != 0:
                    if abs(end - start - length) <= abs(temp_end - temp_start - length):
                        outfile.write(f'{chr}\t{start}\t{end}\tadded_{count_up}_bases\t\t{strand}\n')
                    else:
                        outfile.write(f'{chr}\t{temp_start}\t{temp_end}\tremoved_{count_down}_bases\t\t{strand}\n')
                else:
                    if count_up <= count_down:
                        outfile.write(f'{chr}\t{start}\t{end}\tadded_{count_up}_bases\t\t{strand}\n')
                    else:
                        outfile.write(f'{chr}\t{temp_start}\t{temp_end}\tremoved_{count_down}_bases\t\t{strand}\n')

    # close files
    outfile.close()
    notfile.close()
    openFasta.close()

    print("Finished!")

if __name__ =="__main__":
    args = get_args()
    main(args.inp, args.seq, args.size, args.build, args.outdir, args.thres, args.length, args.up, args.down)
    