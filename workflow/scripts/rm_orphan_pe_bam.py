#!/usr/bin/env python

# source: https://github.com/nf-core/chipseq/blob/master/bin/bampe_rm_orphan.py

# MIT License
#
# Copyright (c) Philip Ewels
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

###############################################################################

#!# own changes and adjustments for snakemake-workflow chipseq are marked with "#!# AVI: " in comment


###############################################################################
###############################################################################
## Created on February 1st 2017 to remove singletons from paired-end BAM file
###############################################################################
###############################################################################

import os
import pysam
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Remove singleton reads from paired-end BAM file i.e if read1 is present in BAM file without read 2 and vice versa.'
Epilog = """Example usage: bampe_rm_orphan.py <BAM_INPUT_FILE> <BAM_OUTPUT_FILE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('BAM_INPUT_FILE', help="Input BAM file sorted by name.")
argParser.add_argument('BAM_OUTPUT_FILE', help="Output BAM file sorted by name.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fr', '--only_fr_pairs', dest="ONLY_FR_PAIRS", help="Only keeps pairs that are in FR orientation on same chromosome.",action='store_true')
args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):

    if not len(path) == 0:
        try:
            #!# AVI: changed because of race conditions if directory exists, original code:  os.makedirs(path)
            os.makedirs(path, exist_ok=True)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def bampe_rm_orphan(BAMIn,BAMOut,onlyFRPairs=False):

    ## SETUP DIRECTORY/FILE STRUCTURE
    OutDir = os.path.dirname(BAMOut)
    makedir(OutDir)

    ## COUNT VARIABLES
    totalReads = 0; totalOutputPairs = 0; totalSingletons = 0; totalImproperPairs = 0

    ## ITERATE THROUGH BAM FILE
    EOF = 0
    SAMFin = pysam.AlignmentFile(BAMIn,"rb")  #!# AVI: changed to new API from pysam.Samfile
    SAMFout = pysam.AlignmentFile(BAMOut, "wb",header=SAMFin.header)   #!# AVI: changed to new API from pysam.Samfile
    currRead = next(SAMFin)     #!# AVI: adapted for the use of the iterator, original code: currRead = SAMFin.next()

    for read in SAMFin.fetch(until_eof=True): #!# AVI: added .fetch() to explicitly use new API
        totalReads += 1
        if currRead.qname == read.qname:
            pair1 = currRead; pair2 = read

            ## FILTER FOR READS ON SAME CHROMOSOME IN FR ORIENTATION
            if onlyFRPairs:
                if pair1.tid == pair2.tid:

                    ## READ1 FORWARD AND READ2 REVERSE STRAND
                    if not pair1.is_reverse and pair2.is_reverse:
                        if pair1.reference_start <= pair2.reference_start:
                            totalOutputPairs += 1
                            SAMFout.write(pair1)
                            SAMFout.write(pair2)
                        else:
                            totalImproperPairs += 1

                    ## READ1 REVERSE AND READ2 FORWARD STRAND
                    elif pair1.is_reverse and not pair2.is_reverse:
                        if pair2.reference_start <= pair1.reference_start:
                            totalOutputPairs += 1
                            SAMFout.write(pair1)
                            SAMFout.write(pair2)
                        else:
                            totalImproperPairs += 1

                    else:
                        totalImproperPairs += 1
                else:
                    totalImproperPairs += 1
            else:
                totalOutputPairs += 1
                SAMFout.write(pair1)
                SAMFout.write(pair2)

            ## RESET COUNTER
            try:
                totalReads += 1
                currRead = next(SAMFin)   #!# AVI: adapted for the use of the iterator, original code: currRead = SAMFin.next()
            except:
                StopIteration
                EOF = 1

        ## READS WHERE ONLY ONE OF A PAIR IS IN FILE
        else:
            totalSingletons += 1
            pair1 = currRead
            currRead = read

    if not EOF:
        totalReads += 1
        totalSingletons += 1
        pair1 = currRead

    ## CLOSE ALL FILE HANDLES
    SAMFin.close()
    SAMFout.close()

    LogFile = os.path.join(OutDir,'%s_bampe_rm_orphan.log' % (os.path.basename(BAMOut[:-4])))
    SamLogFile = open(LogFile,'w')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('FILES/DIRECTORIES')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Input File: ' + BAMIn + '\n')
    SamLogFile.write('Output File: ' + BAMOut + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('OVERALL COUNTS')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Total Input Reads = ' + str(totalReads) + '\n')
    SamLogFile.write('Total Output Pairs = ' + str(totalOutputPairs) + '\n')
    SamLogFile.write('Total Singletons Excluded = ' + str(totalSingletons) + '\n')
    SamLogFile.write('Total Improper Pairs Excluded = ' + str(totalImproperPairs) + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

bampe_rm_orphan(BAMIn=args.BAM_INPUT_FILE,BAMOut=args.BAM_OUTPUT_FILE,onlyFRPairs=args.ONLY_FR_PAIRS)

############################################
############################################
############################################
############################################
