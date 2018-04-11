'''
November 24, 2017, cjyoon@kaist.ac.kr
wrapper script to process fastq files into a  bam file, adapted from ayh's dnaautoaligner.py
currently only supports BWA alignment for WGS, but can be adapted to RNA-seq alignment
After aligning into a bam, uses Picard's mark duplicate, GATK's indel realigner and base recalibration for final processed bam.

Jan 4, 2018 added cleanup module which included option to keep intermediate files per KISTI's request
'''

import os
import subprocess
import sys
import yaml
import argparse
import shlex
import re

def configPath(configurationPath):
    '''configures the path to softwares and input files from a yaml file'''
    with open(configurationPath, 'r') as f:
        pathMaps = yaml.safe_load(f)

    try:
        #print(pathMaps)
        JAVAPATH = pathMaps['JAVAPATH']
        PICARDPATH = pathMaps['PICARDPATH']
        GATKPATH = pathMaps['GATKPATH']
        SAMTOOLSPATH = pathMaps['SAMTOOLSPATH']
        BWAPATH = pathMaps['BWAPATH']
        knownIndelPath = pathMaps['knownIndelPath']
        dbSnpPath = pathMaps['dbSnpPath']
        referencePath = pathMaps['referencePath']
    except KeyError:
        print('one of the required input path not specified. Exiting...')

    return JAVAPATH, PICARDPATH, GATKPATH, SAMTOOLSPATH, BWAPATH,  knownIndelPath, dbSnpPath, referencePath

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', default=os.path.dirname(os.path.realpath(__file__)) + '/config.yml', help='yaml file that can configure absolute paths to executable and databases')
    parser.add_argument('-m', '--memory', default=8, help='memory to allocate for Picard Java operation')
    parser.add_argument('-s', '--sampleID', required=True, help='sample name to be added to the header with @RG SM:')
    parser.add_argument('-t', '--thread', default=4, help='number of threads to be used for samtools processes')
    parser.add_argument('-n', '--ncore', default=4, help='number of multicore to utilize')
    parser.add_argument('-f1', '--fastq1', required=True, help='FASTQ Read1 path')
    parser.add_argument('-f2', '--fastq2', required=True, help='FASTQ Read2 path')
    parser.add_argument('-o', '--outputDIR', required=False, default=os.getcwd(), help='path to which aligned bams will be written')
    #parser.add_argument('-r', '--referencePath', required=True, help='reference FASTA path')
    parser.add_argument('-d', '--dryrun', required=False, default=False, help='print out the command lines, but do not actually run those', type=bool)
    parser.add_argument('--clean', required=False, default=1, choices=[0, 1], type=int, help='If set to 0, then keep all intermediate files, if set to 1, then clean up all intermediate files')
    args = parser.parse_args()

    return args.config, args.fastq1, args.fastq2, args.sampleID, args.memory, args.thread, args.ncore, args.outputDIR, args.dryrun, args.clean

def execute(commandline, dryrun=False):
    ''' Run the command line string as shell script. 
    if dryrun==True, then just print out the command without actually running the command'''
    
    if dryrun:
        print(commandline)
    else:
        execute = subprocess.Popen(shlex.split(commandline))
        execute.wait()

    return 0

def indexBam(SAMTOOLSPATH, bamPath, dryrun=False):
    '''indexes the given bam'''
    indexBamCMD = f'{SAMTOOLSPATH} index -@ 4 {bamPath}'
    execute(indexBamCMD, dryrun)

    return bamPath

def align_sort(fq1, fq2, sampleName, referencePath, BWAPATH, SAMTOOLSPATH, outputDIR, nthread=4, dryrun=False, clean=1):
    ''' performs BWA mem alignment and sorting'''
    alignedSam = f'{outputDIR}/{sampleName}.sam'
    alignedBam = f'{outputDIR}/{sampleName}.bam'
    alignedSortedBam = f'{outputDIR}/{sampleName}.sorted.bam'

    alignCMD = f'{BWAPATH} mem -v 1 -t {nthread} -R "@RG\\tID:{sampleName}\\tSM:{sampleName}\\tPL:ILLUMINA" {referencePath} {fq1} {fq2}'
    if os.path.isfile(alignedBam):
        pass
    else:
        # align to sam file
        if dryrun:
            print(alignCMD)
        else:
            with open(alignedSam, 'w') as f:
                alignExecute = subprocess.Popen(shlex.split(alignCMD), stdout=f)
                alignExecute.wait()

            print('### Done aligning to SAM')
        # convert sam to bam 
        sam2bamCMD = f'{SAMTOOLSPATH} view -bS -O BAM -o {alignedBam} {alignedSam}'
        execute(sam2bamCMD, dryrun)
        if clean==1:
            cleanup([alignedSam], dryrun)

        print('### Done converting SAM to BAM')

	# sort bam 
        bamsortCMD = f'{SAMTOOLSPATH} sort -T {outputDIR}/{sampleName}.sorting -@ {nthread} -O bam -o {alignedSortedBam} {alignedBam}'
        execute(bamsortCMD, dryrun)
        if clean==1:
            cleanup([alignedBam], dryrun)
        print('### Done sorting BAM')

        # index aligned and sorted bam     
        indexBam(SAMTOOLSPATH, alignedSortedBam, dryrun)
        print('### Done indexing aligned and sorted BAM')

    return os.path.abspath(alignedSortedBam)

def markduplicate(JAVAPATH, PICARDPATH, bamPath, dryrun, clean=1):
    '''mark duplicate of the input bam'''
    markedBam =  re.sub(string=bamPath, pattern=r'.bam$', repl='.md.bam')
    # input index bai
    inputBai = f'{bamPath}.bai'

    metricTxt =  re.sub(string=bamPath, pattern=r'.bam$', repl='.md.txt')
    tmp_dir = os.path.dirname(bamPath) 

    mdCMD = f'{JAVAPATH} -XX:ParallelGCThreads=8 -Xmx8g -jar {PICARDPATH} MarkDuplicates REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={bamPath} O={markedBam} M={metricTxt} VALIDATION_STRINGENCY=LENIENT TMP_DIR={tmp_dir} QUIET=true'

    if os.path.isfile(markedBam):
        pass
    else:
        execute(mdCMD, dryrun)
        if clean==1:
            cleanup([bamPath, inputBai], dryrun)
        print('### Done marking duplicates')

    return os.path.abspath(markedBam)

def realign(JAVAPATH, GATKPATH, referencePath, bamPath, knownIndelPath, dryrun, clean=1):
    '''creates indel realinger target from known indels -> realigns the provided bam '''

    realignTargetInterval =  re.sub(string=bamPath, pattern=r'.bam$', repl='.intervals')
    realignedBam = re.sub(string=bamPath, pattern=r'.bam$', repl='.indel.bam')
    realignedBai = re.sub(string=bamPath, pattern=r'.bam$', repl='.indel.bai')

    # input index bai
    inputBai = f'{bamPath}.bai'

    if os.path.isfile(realignedBam):
        pass
    else:

        realignIntervalCMD = f'{JAVAPATH} -Xmx4g -jar {GATKPATH} -T RealignerTargetCreator -R {referencePath} -I {bamPath} --known {knownIndelPath} -o {realignTargetInterval}'
        execute(realignIntervalCMD, dryrun)

        realignBamCMD = f'{JAVAPATH} -Xmx4g -jar {GATKPATH} -T IndelRealigner -R {referencePath} -I {bamPath} -targetIntervals {realignTargetInterval} -o {realignedBam}'
        execute(realignBamCMD, dryrun)
        if clean==1:
            cleanup([bamPath, inputBai], dryrun)

        print('### Done realigning around indels')

    return os.path.abspath(realignedBam)

def baserecalibrator(JAVAPATH, GATKPATH, referencePath, bamPath, dbSnpPath, knownIndelPath, dryrun, clean=1):
    '''recalibrates base quality'''
    
    recalibrateTable = re.sub(string=bamPath, pattern=r'.bam$', repl='.table')
    recalibratedBam = re.sub(string=bamPath, pattern=r'.bam$', repl='.br.bam')

    # input index bai
    inputBai = re.sub(string=bamPath, pattern=r'.bam$', repl='.bai')

    if os.path.isfile(recalibratedBam):
        pass
    else:
        baserecalibrateTableCMD = f'{JAVAPATH} -Xmx4g -jar {GATKPATH} -T BaseRecalibrator -R {referencePath} -I {bamPath} -knownSites {dbSnpPath} --knownSites {knownIndelPath} -o {recalibrateTable}'
        execute(baserecalibrateTableCMD, dryrun)

        baserecalibrateCMD = f'{JAVAPATH} -Xmx4g -jar {GATKPATH} -T PrintReads -R {referencePath} -I {bamPath} -BQSR {recalibrateTable} -o {recalibratedBam} -nct 8'
        execute(baserecalibrateCMD, dryrun)
        if clean==1:
            cleanup([bamPath, inputBai], dryrun)
        print('### Done recalibrating base quality')

    return os.path.abspath(recalibratedBam)

def cleanup(intermediateFileList, dryrun):
    '''removes all the intermediate files that are not necesssary'''
    for f in intermediateFileList:
        cleanupCMD = 'rm -rf ' + f
        execute(cleanupCMD, dryrun)

    return 0

def main():
    config, fastq1, fastq2, sampleID, memory, thread, ncore, outputDIR, dryrun, clean = argument_parser()

    # Configure executale and database paths
    JAVAPATH, PICARDPATH, GATKPATH, SAMTOOLSPATH, BWAPATH, knownIndelPath, dbSnpPath, referencePath = configPath(config)

    # BWA MEM align and sort
    alignedBam = align_sort(fastq1, fastq2, sampleID, referencePath, BWAPATH, SAMTOOLSPATH, outputDIR, thread, dryrun, clean)
    # mark duplicate
    markedBam = markduplicate(JAVAPATH, PICARDPATH, alignedBam, dryrun, clean)
    indexBam(SAMTOOLSPATH, markedBam, dryrun)

    # indel realign
    realignedBam = realign(JAVAPATH, GATKPATH, referencePath, markedBam, knownIndelPath, dryrun, clean)

    # base recalibrate
    recalibratedBam = baserecalibrator(JAVAPATH, GATKPATH, referencePath, realignedBam, dbSnpPath, knownIndelPath, dryrun, clean)

    if dryrun:
        print('### This was a dry run')
    else:
        print(f'### Analysis ready bam is located {os.path.abspath(recalibratedBam)}')

    return 0

if __name__ == '__main__':
    main()

