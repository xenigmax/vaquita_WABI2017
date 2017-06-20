#!/usr/bin/python
'''
RUN

Jongkyu Kim@MPI-Molgen/FU-Berlin
xenigmax@gmail.com
'''

import os
import sys
import argparse
import time
import subprocess

sys.path.append(os.getcwd())
import config as c 

def v_seqan(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/v_seqan#" + ("#".join(param[3:])).replace("/","_")
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    logfile_name = outdir + "/log.txt"
    outfile_name = outdir + "/out.vcf"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["v_seqan"]) # bin
    for option in param[3:] : # options
        if option[0] == "-" and option[1] != "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
           args.append(option)
    args.append(bamdir + "/out.sorted.bam")

    print("[V_SEQAN] START : %s outdir=%s, logfile=%s" % (" ".join(args), outdir, logfile_name))
    try :
        checkAndMake(outdir)
        outfile = open(outfile_name, "w") 
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=outfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[V_SEQAN] ERROR")
        #print e.output
        return
    print("[V] END")


def v(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/v#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    logfile_name = outdir + "/log.txt"

    # args
    bUseRealign = False
    args = []
    args.append(c.CONFIG["variant-callers"]["v"]) # bin
    args.append(bamdir + "/stellar-gff.gff")
    for option in param[3:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            if option == "realign" :
                bUseRealign = True
            else :
                args.append(option)
    args.append(outdir)

    # realign
    if bUseRealign : 
        inbam = bamdir + "/out.sorted.bam" 
        outsam = outdir + "/out.unmapped.sam"
        outfa = outsam.replace(".sam", ".fa")
        extractUnmappedReads(inbam, outsam)
        unmappedSam2Fa(outsam, outfa)

        outgff = outsam.replace(".sam", ".gff")
        _param = []
        _param.append( param[0] )
        _param.append( param[1] )
        _param.append( "-l:20" )
        _param.append( "-k:16" )
        _param.append( "-e:0.03" )
        stellar(outfa, outgff, _param)
 
    exit(1)
    print("[V] START : %s outdir=%s, logfile=%s" % (" ".join(args), outdir, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[V] ERROR")
        #print e.output
        return
    print("[V] END")

def gustaf(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/gustaf#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    outvcf_name = outdir + "/out.vcf"  
    outgff_name = outdir + "/out.gff"  
    logfile_name = outdir + "/log.txt"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")    
    args.append(c.CONFIG["variant-callers"]["gustaf"]) # bin
    args.append(c.CONFIG["ref-dir"] + "/" + ref) # ref.fa

    reads = reads.split(",")
    if len(reads) > 1 : # gustaf can't load multiple files.
        inFiles = []
        for read in reads :
            read = read.replace(".fq",".fa")
            inFiles.append(c.CONFIG["read-dir"] + "/" + read) 
        outFile = outdir + "/merged.fa"
        print inFiles, outFile
        mergeFiles(inFiles, outFile)
        args.append(outFile)
    else :
        args.append(c.CONFIG["read-dir"] + "/" + reads[0]) 

    args.append("-m")
    args.append(bamdir + "/stellar-gff.gff")
    for option in param[3:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)
    args.append("-gff")
    args.append(outgff_name)
    args.append("-vcf")
    args.append(outvcf_name)

    print("[GUSTAF] START : %s outfile=%s, logfile=%s" % (" ".join(args), outvcf_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[GUSTAF] ERROR")
        print e.output
        return
    print("[GUSTAF] END")

def pindel(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/pindel#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    outfile_name = outdir + "/out.txt"  
    logfile_name = outdir + "/log.txt"

    # make inputfile
    inputfile_name = outdir + "/in.txt"
    infile = open(inputfile_name,"w")
    infile.write("%s %s %s\n" % (bamdir + "/out.sorted.bam", param[3], "out")) # param[3] : insertion size
    infile.close()

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["pindel"]) # bin
    args.append("-f")
    args.append(c.CONFIG["ref-dir"] + "/" + ref) # ref.fa
    args.append("-i")
    args.append(inputfile_name) # ref.fa
    args.append("-o")
    args.append(outfile_name)
    for option in param[4:]: # options (param[3] : insertion size)
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)

    print("[PINDEL] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[PINDEL] ERROR")
        #print e.output
        return
    print("[PINDEL] END")

def delly(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/delly#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    outfile_name = outdir + "/out.bcf"  
    logfile_name = outdir + "/log.txt"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["delly"]) # bin
    for option in param[3:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)
    args.append("-g")
    args.append(c.CONFIG["ref-dir"] + "/" + ref) # ref.fa
    args.append("-o")
    args.append(outfile_name)
    args.append(bamdir + "/out.sorted.bam")

    print("[DELLY] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[DELLY] ERROR")
        #print e.output
        return
    print("[DELLY] END")
    bcf2vcf(outfile_name)

def lumpy(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

    outdir = bamdir + "/lumpy#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    outfile_name = outdir + "/out.vcf"  
    logfile_name = outdir + "/log.txt"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["lumpy"]) # bin
    for option in param[3:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)
    args.append("-B")
    args.append(bamdir + "/out.sorted.LUMPY.bam")
    #args.append(bamdir + "/out.sorted.bam")
    args.append("-T")
    args.append(bamdir + "/temp_LUMPY")
    args.append("-o")
    args.append(outfile_name)
    args.append("-k")

    print("[LUMPY] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[LUMPY] ERROR")
        #print e.output
        return
    print("[LUMPY] END")


def gasvpro(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = param[2] # .sorted.bam

        #/srv/public/kjk/noname/bin/gasv/bin/GASVPro.sh ../out.sorted.hq.GASV.bam ../out.name_sorted.lq.GASV.bam ./temp_GASV


    outdir = bamdir + "/gasvpro#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    logfile_name = outdir + "/log.txt"
    outfile_name ="program's default"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["gasvpro"]) # bin
    args.append("../out.sorted.hq.GASV.bam")
    args.append("../out.sorted.lq.GASV.bam")
    args.append("./temp_GASV")
    args.append(outdir)

    print("[GASVPro] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[GASVPro] ERROR")
        #print e.output
        return
    print("[GASVPro] END")

def crest(param, getoutdir=False) :
    ref = param[0] # ref.fa
    reads = param[1]
    bamdir = "/srv/public/kjk/noname/analysis/" + param[2] # .sorted.bam, but full path

    outdir = bamdir + "/crest#" + "#".join(param[3:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    logfile_name = outdir + "/log.txt"
    outfile_name ="program's default"
    genome_fa_file = "/srv/public/kjk/noname/data/refs/" + ref;
    genome_2bit_file = genome_fa_file.replace(".fa",".2bit")

    if "22" in ref :
        blat_port = "6622"
    else :
        blat_port = "6619"

    # args
    args = []
    args.append("/usr/bin/time")
    args.append("-v")
    args.append(c.CONFIG["variant-callers"]["crest"]) # bin
    args.append(blat_port)
    args.append(bamdir + "/out.sorted.bam")
    args.append(genome_fa_file)
    args.append(genome_2bit_file)
    args.append(outdir)
    print args

    print("[crest] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[crest] ERROR")
        #print e.output
        return
    print("[crest] END")

def bcf2vcf(filename) :
    print("[BCF2VCF] START")
    try :
        outfile = open(filename.replace(".bcf",".vcf"), "w") 
        args = []
        args.append(c.CONFIG["tools"]["bcftools"])
        args.append("view")
        args.append(filename)
        subprocess.check_call(args, stdout=outfile)
        outfile.close()
    except subprocess.CalledProcessError as e: 
        print("[BCF2VCF] ERROR")
        return
    print("[BCF2VCF] END")


def stellar(infile_name, outfile_name, param) :
    print param

    ref = param[0] # reference.fa

    args = []
    args.append(c.CONFIG["read-mappers"]["stellar"]) 
    args.append(c.CONFIG["ref-dir"] + "/" + ref) # ref.fa
    args.append(infile_name) 
    for option in param[2:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)
    args.append("-o") 
    args.append(outfile_name)
    logfile_name = outfile_name + ".log"
    print outfile_name
    print logfile_name

    print("[STELLAR] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=logfile, stderr=logfile)
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[STELLAR] ERROR")
        #print e.output
        return
    print("[STELLAR] END")


def bwa(param, getoutdir=False) :
    ref = param[0] # reference.fa
    reads = param[1]

    outdir = c.CONFIG["result-dir"] + "/" + ref + "/" + reads.replace(",","#") + "/bwa#" + "#".join(param[2:])
    if getoutdir :
        return outdir

    checkAndMake(outdir)
    outfile_name = outdir + "/out.sam"  
    logfile_name = outdir + "/log.txt"

    # args
    args = []
    args.append(c.CONFIG["read-mappers"]["bwa"]) # bin
    for option in param[2:] : # options
        if option[0] == "-" :
            s = option.split(":")
            args.append(s[0])
            args.append(s[1])
        else :
            args.append(option)
    args.append(c.CONFIG["ref-dir"] + "/" + ref) # ref.fa
    for read in reads.split(",") :
        args.append(c.CONFIG["read-dir"] + "/" + read) # reads

    print("[BWA] START : %s outfile=%s, logfile=%s" % (" ".join(args), outfile_name, logfile_name))
    try :
        checkAndMake(outdir)
        outfile = open(outfile_name, "w") 
        logfile = open(logfile_name, "w") 
        subprocess.check_call(args, stdout=outfile, stderr=logfile)
        outfile.close()
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[BWA] ERROR")
        #print e.output
        return
    print("[BWA] END")

    # sam2bam
    sam2bam(outfile_name)
    sortBam(outfile_name.replace(".sam",".bam"))
    indexBam(outfile_name.replace(".sam",".sorted.bam"))

    # sam2stellargff
    #sam2stellargff(outfile_name)

    return outdir

# extract unmapped reads(filename)
def extractUnmappedReads(inbam, outsam) :
    print("[EXTRACT UNMAPPED READS] START")
    try :
        args = []
        args.append(c.CONFIG["tools"]["samtools"])
        args.append("view")
        args.append("-f")
        args.append("4") # flag for unmappared reads
        args.append(inbam)
        outfile = open(outsam, "w") 
        subprocess.check_call(args, stdout=outfile)
        outfile.close()
    except subprocess.CalledProcessError as e: 
        print("[EXTRACT UNMAPPED READS] ERROR")
        return
    print("[EXTRACT UNMAPPED READS] END")

def unmappedSam2Fa(insam, outfa) :
    infile = open(insam, "r")
    outfile = open(outfa, "w")
    
    for line in infile :
        s = line.strip().split()
        readname = s[0]
        flag = int(s[1])
        seq = s[9]

        if (flag & 0x40) > 0  : # first in pair
            readname += "/1"
        if (flag & 0x80) > 0 : # second in pair
            readname += "/2"
        outfile.write(">%s\n%s\n" % (readname, seq))
    infile.close()
    outfile.close()

# .sam to .bam
def sam2bam(filename) :
    print("[SAM2BAM] START")
    try :
        args = []
        args.append(c.CONFIG["tools"]["samtools"])
        args.append("view")
        args.append("-b")
        args.append("-S")
        args.append(filename)
        args.append("-o")
        args.append(filename.replace(".sam",".bam")) # .sam => .bam
        subprocess.check_call(args)
    except subprocess.CalledProcessError as e: 
        print("[SAM2BAM] ERROR")
        return
    print("[SAM2BAM] END")

# .bam to .sorted.bam
def sortBam(filename) :
    print("[SORT-BAM] START")
    try :
        args = []
        args.append(c.CONFIG["tools"]["samtools"])
        args.append("sort")
        args.append(filename)
        args.append("-o")
        args.append(filename.replace(".bam",".sorted.bam")) # .bam => .sorted.bam
        subprocess.check_call(args)
    except subprocess.CalledProcessError as e: 
        print("[SORT-BAM] ERROR")
        return
    print("[SORT-BAM] END")

# make indeX
def indexBam(filename) :
    print("[INDEX-BAM] START")
    try :
        args = []
        args.append(c.CONFIG["tools"]["samtools"])
        args.append("index")
        args.append(filename) # .bam => .bam.bai
        subprocess.check_call(args)
    except subprocess.CalledProcessError as e: 
        print("[INDEX-BAM] ERROR")
        return
    print("[INDEX-BAM] END")


# make stellar-compatiable .gff (for Gustaf..)
def sam2stellargff(filename) :
    print("[SAM2STELLAR-GFF] START")
    try :
        dirname = os.path.dirname(filename)
        outfile = open(dirname + "/stellar-gff.gff", "w") 
        logfile = open(dirname + "/stellar-gff.within-indel", "w") 
        args = []
        args.append("perl")
        args.append(c.CONFIG["tools"]["sam2stellargff"])
        args.append(filename)
        subprocess.check_call(args, stdout=outfile, stderr=logfile)
        outfile.close()
        logfile.close()
    except subprocess.CalledProcessError as e: 
        print("[SAM2STELLAR-GFF] ERROR")
        return
    print("[SAM2STELLAR-GFF] END")

def mergeFiles(inFiles, outFilename) :
    print("[MERGEFILES] START")
    try :
        outfile = open(outFilename, "w") 
        args = []
        args.append("cat")
        for filename in inFiles :
            args.append(filename)
        subprocess.check_call(args, stdout=outfile)
        outfile.close()
    except subprocess.CalledProcessError as e: 
        print("[MERGEFILES] ERROR")
        return
    print("[MERGEFILES] END")

#########################################################################
class Logger(object) :
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def __del__(self):
        self.log.close()

    def write(self, message):
        if message.strip() != "" :
            timetag = time.strftime("[%Y%m%d %H:%M:%S] ") 
            _msg = timetag + message + "\n"
            self.terminal.write(_msg)
            self.log.write(_msg)

def checkAndMake(dirname) :
    if not os.path.exists(dirname) :
        try :
            os.makedirs( dirname )
        except Exception:
            print("failed to make [%s]", dirname)

def checkDirs():
    checkAndMake(c.CONFIG["result-dir"])
    checkAndMake(c.CONFIG["log-dir"])

def run(func, param1, param2) :
    return eval(func)(param1, param2)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--reference", help="reference file in FASTA format", required=True)
    parser.add_argument("-r","--reads", help="read files(delimeter : ,) in FASTA or FASTQ format", required=True)

    rm_names = "|".join([key for key in c.CONFIG["read-mappers"]])
    parser.add_argument("-m","--read-mapper", help="[%s],option1,option2,.." % rm_names, required=True)
    vc_names = "|".join([key for key in c.CONFIG["variant-callers"]])
    parser.add_argument("-v","--variant-caller", help="[%s],option1,option2,.." % vc_names, required=False)
    parser.add_argument("--skip-read-mapping", help="skip this step", action='store_true', required=False)
    parser.add_argument("--skip-variant-calling", help="skip this step", action='store_true', required=False)
    args = parser.parse_args()

    # check dirs
    checkDirs()

    # open logger
    sys.stdout = Logger(c.CONFIG["log-dir"] + "/" + time.strftime("%Y%m%d.txt"))
    print("[START] " + " ".join(sys.argv))

    # run align
    rm = args.read_mapper.split(",")
    rm_name = rm[0]
    rm_options = [args.reference] + [args.reads] + rm[1:]
    if args.skip_read_mapping :
        print("[FUNC] SKIPPED: %s" % args.read_mapper)
        outdir = run(rm_name, rm_options, True )
    else :
        if (rm_name in dir()) == False :
            print("[ERROR] Can't find the function call : %s" % rm_name)
            exit(1)
        else :
            print("[FUNC] START: %s" % args.read_mapper)
            outdir = run(rm_name, rm_options, False)
            print("[FUNC] END: %s" % args.read_mapper)

    # run vc
    if args.skip_variant_calling :
        if args.variant_caller :
            print("[FUNC] SKIPPED: %s" % args.variant_caller)
        else :
            print("[FUNC] SKIPPED: VARIANT CALLING")
    else :
        vc = args.variant_caller.split(",")
        vc_name = vc[0]
        vc_options = [args.reference] + [args.reads] + [outdir] +  vc[1:]
        if (vc_name in dir()) == False :
            print("[ERROR] Can't find the function call : %s" % rm_name)
            exit(1)
        else :
            print("[FUNC] START: %s" % args.variant_caller)
            run(vc_name, vc_options, False)
            print("[FUNC] END: %s" % args.variant_caller)

    print("[END]")
