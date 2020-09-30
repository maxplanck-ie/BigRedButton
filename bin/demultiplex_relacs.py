#!/usr/bin/env python
import sys
import os
import argparse
import glob
from multiprocessing import Pool
import subprocess
import gzip
import matplotlib.pyplot as plt
import numpy as np

try:
    import editdistance as ed
except:
    sys.exit("You must install 'editdistance' first")


def parseArgs(args=None):
    parser = argparse.ArgumentParser(description="Process RELACS reads by moving a barcode into the header and writing output to per-barcode files.")
    parser.add_argument('-p', '--numThreads', help='Number of threads to use. Note that there will only ever be a single thread used per illumina sample, so there is no reason so specify more threads than samples. Note also that the effective number of threads used will always be higher than this, since compression/decompression for each file occurs in a separate thread.', type=int, default=1)
    parser.add_argument('-b', '--buffer', help='Number of bases of "buffer" to ignore after the barcode. (default: 1)', type=int, default=1)
    parser.add_argument('--umiLength', type=int, default=0, help="If present, the number of bases to be used as the UMI. These proceed the cell barcode.")
    parser.add_argument('sampleTable', help="""A tab-separated table with three columns: Illumina sample	barcode	sample_name

An example is:

Sample_Lib_C838_11	ACTACT	Sample1
Sample_Lib_C838_11	TGACTG	Sample2
Sample_Lib_C838_11	default	Sample3

The "Illumina sample" column corresponds to a directory with fastq files, as produced by bcl2fastq2.

Any observed barcode that does not match a barcode listed in the file (with a possible mismatch of 1) will be written to the sample specified by 'default'. Note that all barcodes must be of the same length.

This script must be run in directory containing subdirectories each having fastq files. This is the default structure output by bcl2fastq2 and our internal demultiplexing pipeline.
""")
    parser.add_argument('output', metavar='output_basename', help="Base output directory, which must exist. As an example, if a given PE read has the barcode ACTACT, then it will be written to /output_basename/Sample_Lib_C838_11/Sample_R1.fastq.gz and /output_basename/Sample_Lib_C838_11/Sample_R2.fastq.gz")
    parser.add_argument('--version', action='version', version="%(prog)s 1.0")
    args = parser.parse_args(args)

    if args.sampleTable is None or args.output is None:
        parser.print_help()

    return args


def readSampleTable(sampleTable):
    """
    Read a sample table into a dict (keys are barcodes).
    Return the resulting dict and the barcode length.
    """
    f = open(sampleTable)
    d = dict()
    bcLen = 0
    for line in f:
        _ = line.strip().split("\t")
        if len(_) < 3:
            continue
        if _[0] not in d:
            d[_[0]] = dict()
        d[_[0]][_[1]] = _[2]
        if _[1] != 'default' and len(_[1]) > bcLen:
            bcLen = len(_[1])
    f.close()

    return (d, bcLen)


def matchSample(sequence, sequence2, oDict, bcLen, umiLength):
    """
    Match a barcode against the sample sheet with a possible edit distance of 1.

    Returns a tuple of the list of file pointers and True/False (whether the "default" output is being used)
    """
    bc = sequence[umiLength:bcLen + umiLength]
    bc2 = sequence2[umiLength:bcLen + umiLength] if sequence2 else None  # For whatever reason, padding isn't used
    if bc in oDict:
        if bc2 and ed.eval(bc, bc2) < 2:
            return (bc, True)
        if not bc2:
            return (bc, True)

    # Look for a 1 base mismatch
    for k, v in oDict.items():
        if ed.eval(k, bc) == 1:
            if not bc2:
                return (k, True)
            if bc2 and ed.eval(bc, bc2) < 2:
                return (k, True)

    return ("default", False)


def writeRead(lineList, of, bc, bcLen, args, doTrim=True):
    """
    Trim the read as specified and write to the appropriate file.

    Return the modified read name (for read #2, if needed)
    """
    rname = lineList[0]
    if doTrim:
        UMI = ""
        if args.umiLength > 0:
            UMI = lineList[1][:args.umiLength]

        # Trim off the barcode
        lineList[1] = lineList[1][bcLen + args.buffer + args.umiLength:]
        lineList[3] = lineList[3][bcLen + args.buffer + args.umiLength:]

        # Fix the read name
        rname = rname.split()
        if args.umiLength > 0:
            rname[0] = "{}_{}_{}".format(rname[0], bc, UMI)
        else:
            rname[0] = "{}_{}".format(rname[0], bc)
        rname = " ".join(rname)

    if rname[-1] != '\n':
        rname += '\n'

    of.write(rname.encode())
    of.write(lineList[1].encode())
    of.write(lineList[2].encode())
    of.write(lineList[3].encode())

    return rname, bc


def writeRead2(lineList, of, bcLen, args, doTrim=True):
    # Fix the read name so it's read #2 rather than #1
    rname = lineList[0]
    rname = rname.split()
    rname[1] = "2{}".format(rname[1][1:])
    rname = " ".join(rname)

    if rname[-1] != '\n':
        rname += '\n'

    if doTrim:
        # trim off the barcode
        lineList[1] = lineList[1][bcLen + args.buffer + args.umiLength:]
        lineList[3] = lineList[3][bcLen + args.buffer + args.umiLength:]

    of.write(rname.encode())
    of.write(lineList[1].encode())
    of.write(lineList[2].encode())
    of.write(lineList[3].encode())

def writePaired(read1, read2, of, bc, bcLen, args, doTrim=True):
    """
    """
    rname = read1[0]
    if doTrim:
        UMI = ""
        if args.umiLength > 0:
            UMI = read1[1][:args.umiLength]
            UMI += read2[1][:args.umiLength]

        # Trim off the barcode
        read1[1] = read1[1][bcLen + args.buffer + args.umiLength:]
        read1[3] = read1[3][bcLen + args.buffer + args.umiLength:]

        read2[1] = read2[1][bcLen + args.buffer + args.umiLength:]
        read2[3] = read2[3][bcLen + args.buffer + args.umiLength:]

        # Fix the read name
        rname = rname.split()
        if args.umiLength > 0:
            rname[0] = "{}_{}_{}".format(rname[0], bc, UMI)
        else:
            rname[0] = "{}_{}".format(rname[0], bc)
        rname = " ".join(rname)

    if rname[-1] != '\n':
        rname += '\n'
    
    of[0].write(rname.encode())
    of[0].write(read1[1].encode())
    of[0].write(read1[2].encode())
    of[0].write(read1[3].encode())

    of[1].write(rname.encode())
    of[1].write(read2[1].encode())
    of[1].write(read2[2].encode())
    of[1].write(read2[3].encode())

    return bc


def processPaired(args, sDict, bcLen, read1, read2, bc_dict):
    f1_ = subprocess.Popen("gunzip -c {}".format(read1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    f2_ = subprocess.Popen("gunzip -c {}".format(read2), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    f1 = f1_.stdout
    f2 = f2_.stdout
    false_bc = 0
    for line1_1 in f1:
        line1_1 = line1_1.decode("ascii")
        line1_2 = f1.readline().decode("ascii")
        line1_3 = f1.readline().decode("ascii")
        line1_4 = f1.readline().decode("ascii")
        line2_1 = f2.readline().decode("ascii")
        line2_2 = f2.readline().decode("ascii")
        line2_3 = f2.readline().decode("ascii")
        line2_4 = f2.readline().decode("ascii")
        (bc, isDefault) = matchSample(line1_2, line2_2, sDict, bcLen, args.umiLength)

        relacs_bc = writePaired([line1_1, line1_2, line1_3, line1_4], [line1_2,line2_2, line2_3, line2_4], sDict[bc], bc, bcLen, args, isDefault)

        if isDefault is True:
           if relacs_bc not in bc_dict.keys():
              bc_dict[bc] = 1
           else:
              bc_dict[bc] += 1
        else: 
              false_bc += 1
    print(bc_dict, false_bc)
    plot_bc_occurance(read1, bc_dict, false_bc, args.output)
    
#   rname, relacs_bc = writeRead([line1_1, line1_2, line1_3, line1_4], sDict[bc][0], bc, bcLen, args, isDefault)
#   writeRead2([rname , line2_2, line2_3, line2_4], sDict[bc][1], bcLen, args, isDefault)

    f1.close()
    f2.close()


def processSingle(args, sDict, bcLen, read1):
    f1_ = subprocess.Popen("gunzip -c {}".format(read1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    f1 = f1_.stdout

    for line1_1 in f1:
        line1_1 = line1_1.decode("ascii")
        line1_2 = f1.readline().decode("ascii")
        line1_3 = f1.readline().decode("ascii")
        line1_4 = f1.readline().decode("ascii")
        (bc, isDefault) = matchSample(line1_2, None, sDict, bcLen, args.umiLength)
        writeRead([line1_1, line1_2, line1_3, line1_4], sDict[bc], bc, bcLen, args, isDefault)

    f1.close()


def wrapper(foo):
    d, args, sDict, bcLen,bc_dict = foo
    print("Processing library {}".format(d))

    # Make the output directories
    try:
        os.makedirs("{}/{}".format(args.output, d))
    except:
        pass

    # Find the fastq files
    R1 = glob.glob("{}/*_R1.fastq.gz".format(d))
    R2 = None
    if len(R1) > 1:
        print("Warning, there was more than 1 sample found in {}, only {} will be used".format(d, R1[0]))
    R1 = R1[0]
    if os.path.exists("{}_R2.fastq.gz".format(R1[:-12])):
        R2 = "{}_R2.fastq.gz".format(R1[:-12])

    # Open the output files and process
    oDict = dict()
    if R2 is not None:
        for k, v in sDict.items():
            oDict[k] = [subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R1.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin,
                        subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R2.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin]
        if 'default' not in oDict:
            k = 'default'
            v = 'unknown'
            oDict[k] = [subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R1.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin,
                        subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R2.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin]
        processPaired(args, oDict, bcLen, R1, R2, bc_dict)
    else:
        for k, v in ssDict.items():
            oDict[k] = [subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R1.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin]
        if 'default' not in oDict:
            k = 'default'
            v = 'unknown'
            oDict[k] = [subprocess.Popen(['gzip', '-c'], stdout=open('{}/{}/{}_R1.fastq.gz'.format(args.output, d, v), "wb"), stdin=subprocess.PIPE, bufsize=0).stdin]
        processSingle(args, oDict, bcLen, R1)
    return bc_dict

def plot_bc_occurance(R1, bc_dict, false_bc, output_path):
    total_sum = false_bc
    for k,v in bc_dict.items():
        total_sum += v

    percentages = [float(false_bc/total_sum)*100]
    x_ticks = ["false_bc"]
    for k,v in bc_dict.items():
       percentages.append(float(v/total_sum)*100)
       x_ticks.append(str(k))
    
    percentages = np.asarray(percentages)
    bc_mean = np.mean(percentages[1:])
    exp_value = 100/float(len(percentages[1:]))
    bc_std = np.std(percentages[1:] - exp_value)
    print("mean ", bc_mean, bc_std)

    print(percentages)
    print(x_ticks)
    fig,ax = plt.subplots(dpi=300)
    x = np.arange(len(percentages))
    ax.bar(x, percentages)
    ax.set_xticks(x)
    ax.set_xticklabels(x_ticks, rotation='vertical', fontsize = 6)
    exp_value = 100/(len(x)-1)
    ax.axhline(y=exp_value, linestyle="--", linewidth=0.5, color='k')
    xx = [-1]+list(range(len(x)))+[len(x)+1]
    print(xx)
    ax.fill_between(xx, [bc_mean + bc_std]*len(xx), [bc_mean - bc_std]*len(xx), color='dimgrey', alpha=0.2, zorder=3)
    plt.ylabel("% of total reads")
    sample_name = R1.split("_R1")[0]
    fig_path_name = os.path.join(output_path,sample_name+"_fig.png")
    plt.savefig(fig_path_name, pad_inches=0.6, bbox_inches='tight')
def main(args=None):
    args = parseArgs(args)
    bc_dict = dict()
    sDict, bcLen = readSampleTable(args.sampleTable)
    p = Pool(processes=args.numThreads)
    tasks = [(d, args, v, bcLen, bc_dict) for d, v in sDict.items()]
    this_bc_dict = p.map(wrapper, tasks)

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["-h"]
    main(args)
