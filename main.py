import os
from readAligner import readAlign
from sam_file import Header, SAMLine
from genomeIndexing import reverseComplement
import sys
import gzip
from urllib import request


def readFastQ(filename1, filename2, output):

    with open(filename1, 'r') as f1, open(filename2, 'r') as f2:
        while True:

            # Read attributes for first read

            readName = f1.readline()[1:].split('/')[0]
            refName = readName.split('-')[0]
            f2.readline()

            readSeq1 = f1.readline().rstrip()
            f1.readline()
            readQuality1 = f1.readline().rstrip()

            # Read attributes for second read

            readSeq2 = f2.readline().rstrip()
            f2.readline()
            readQuality2 = f2.readline().rstrip()

            if len(readSeq1) == 0 or len(readSeq2) == 0:
                break

            # Getting optimal alignments for both reads (Doing only single alignment for now)
            costMethod = [0, 3, 1]

            alignments1, pos1, cigar1 = rA.align(readSeq1, costMethod)
            alignments2, pos2, cigar2 = rA.align(readSeq2, costMethod)

            # Output Format Generation for First Read

            strandsAligned = [alignments1[1][0], alignments2[1][0]]
            flag_bits = rA.getFLAGbits(strandsReadPair=strandsAligned, nAlignRead1=1, nAlignRead2=1, first=True)

            SAMOutput = sF.SAMLine()
            SAMOutput.editFLAG(flag_bits)

            pos = [pos1, pos2]
            cigars = [cigar1, cigar2]
            SAMOutput.editTLEN(pos, cigars)

            if strandsAligned[0] == 'F':
                SAMOutput.editFields(readName, refName, cigar1, pos1, MAPQ=99, RNEXT='=', PNEXT=pos2,
                                     SEQ=readSeq1, qualities=readQuality1)
            else:
                readQuality1 = ''.join(item[1] for item in sorted(enumerate(readQuality1), reverse=True))
                SAMOutput.editFields(readName, refName, cigar1, pos1, MAPQ=99, RNEXT='=', PNEXT=pos2,
                                     SEQ=reverseComplement(readSeq1), qualities=readQuality1)

            SAMOutput.writeToFile(output)

            # Output Format Generation for Second Read

            strandsAligned = [alignments2[1][0], alignments1[1][0]]
            flag_bits = rA.getFLAGbits(strandsReadPair=strandsAligned, nAlignRead1=1, nAlignRead2=1, first=False)

            SAMOutput = sF.SAMLine()
            SAMOutput.editFLAG(flag_bits)

            pos = [pos2, pos1]
            cigars = [cigar2, cigar1]
            SAMOutput.editTLEN(pos, cigars)

            if strandsAligned[0] == 'F':
                SAMOutput.editFields(readName, refName, cigar2, pos2, MAPQ=99, RNEXT='=', PNEXT=pos1,
                                     SEQ=readSeq2, qualities=readQuality2)
            else:
                readQuality2 = ''.join(item[1] for item in sorted(enumerate(readQuality2), reverse=True))
                SAMOutput.editFields(readName, refName, cigar2, pos2, MAPQ=99, RNEXT='=', PNEXT=pos1,
                                     SEQ=reverseComplement(readSeq2), qualities=readQuality2)
            SAMOutput.writeToFile(output)

def get_download_files(url, coverage):
    file1 = url + 'output_' + str(coverage) + 'xCov1.fq.gz'
    file2 = url + 'output_' + str(coverage) + 'xCov2.fq.gz'
    return [file1, file2]

if __name__ == "__main__":

    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    url = 'http://public.bmi.inf.ethz.ch/teaching/cbm1_2017/project1/data/'
    genome_file = 'genome.chr22.fa'
    genome_path = os.path.join(DATA_DIR, genome_file)

    if not os.path.exists(genome_path):
        request.urlretrieve(url + genome_file, genome_path)

    read_aligner = readAlign()
    read_aligner.createGenomeIndex(filename=genome_path, params=[20, 10], method='kmer')

    coverages = [5, 10, 30]

    for coverage in coverages:

        files_to_download = get_download_files(url=url, coverage=coverage)

        file1 = os.path.join(DATA_DIR, 'output_' + str(i) + 'xCov1.fq')
        file2 = os.path.join(DATA_DIR, 'output_' + str(i) + 'xCov2.fq')
        output = os.path.join(DATA_DIR, 'output_' + str(i) + 'xCov.sam')

        if not os.path.exists(file1 + '.gz'):
            request.urlretrieve(downloadfile1, file1 + '.gz')

        if not os.path.exists(file2 + '.gz'):
            request.urlretrieve(downloadfile2, file2 + '.gz')

        with gzip.open(file1 + '.gz', 'rb') as f1:
            with open(file1, 'wb') as f2:
                f2.write(f1.read())

        with gzip.open(file2 + '.gz', 'rb') as f1:
            with open(file2, 'wb') as f2:
                f2.write(f1.read())

        print("[INFO]: Starting read alignment for coverage: {}".format(str(i)))

        header_file = Header(seq_name=read_align.genome_name, genome_length=read_align.genome_length, 
            sorting_order="unsorted")
        header_file.write(output)

        readFastQ(file1, file2, output)
