#!/usr/bin/python
'''
REPORT

Jongkyu Kim@MPI-Molgen/FU-Berlin
xenigmax@gmail.com
'''

import os
import sys
import argparse
import time
import subprocess
import datetime
import types
import math
import copy

sys.path.append(os.getcwd())
import config as c 

TIME_ELAPSED_STR = "TIME_ELAPSED"
TIME_CPU_STR = "TIME_CPU"
MEMORY_STR = "MEMORY"
DEL_STR = "DEL"
INV_STR = "INV"
INS_STR = "INS"
TRA_STR = "TRA" # translocation
BND_STR = "BND" # translocation in mason
DUP_STR = "DUP"
DUP_IMPR_STR = "DUP_IMPR_STR" # imprecise duplication
ITX_STR = "ITX" # for CREST
CNV_STR = "CNV"
BREAKPOINT_STR = "BP"
SUPP_INFO_STR = "_SUPP_INFO"
VARIATIONS = [DEL_STR, INS_STR, DUP_STR, DUP_IMPR_STR, INV_STR, TRA_STR, CNV_STR, BREAKPOINT_STR]
MAX_SUPP = 5

####################################
# True postivies
####################################
def varsim(param) :
    p_name = param[0]

    p = {}
    #for t in VARIATIONS :
    for t in  [DEL_STR]:
        p[t] = []

    p_file = open(c.CONFIG["p-dir"] + "/" + p_name, "r")
    for line in p_file :
        if line[0] == "#" :
            continue

        s = line.strip().split()

        chrm = s[0]
        pos = int(s[1])
        varID = s[2]
        info = parseVcfInfo(s[7]) 

        if ("SVLEN" in info) == False :
            continue

        svlen = int(info["SVLEN"].split(",")[0])

        if svlen < 0 :
            startPos = pos 
            endPos = pos + abs(svlen) - 1
            p[DEL_STR].append( [chrm, startPos, endPos ] ) # [start, end], 1-based 
    p_file.close()

    return p

def mason(param) :
    p_name = param[0]

    p = {}
    for t in VARIATIONS :
        p[t] = []

    traBuf = {}
    p_file = open(c.CONFIG["p-dir"] + "/" + p_name, "r")
    for line in p_file :
        if line[0] == "#" :
            continue

        s = line.strip().split()

        chrm = s[0]
        pos = int(s[1])
        varID = s[2]
        info = parseVcfInfo(s[7]) 

        if ("SVTYPE" in info) == False :
            info["SVTYPE"] = "SNP"

        if info["SVTYPE"] == DEL_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            p[DEL_STR].append([chrm, startPos, endPos]) # [start, end], 1-based
            p[BREAKPOINT_STR].append([chrm, startPos, endPos, DEL_STR])
        elif info["SVTYPE"] == INS_STR :
            startPos = pos 
            endPos = pos + int(info["SVLEN"]) - 1
            p[INS_STR].append( [chrm, startPos, endPos ] ) # [start, end], 1-based 
        elif info["SVTYPE"] == INV_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            p[INV_STR].append( [chrm, startPos, endPos ] )
            p[BREAKPOINT_STR].append([chrm, startPos, endPos, INV_STR])
        elif info["SVTYPE"] == DUP_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            target = info["TARGETPOS"].split(":")
            targetChrm = target[0]
            targetStartPos = int(target[1] )
            targetEndPos = targetStartPos + abs(int(info["SVLEN"])) - 1
            p[DUP_STR].append( [chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )

            if startPos < targetStartPos :
                p[DUP_IMPR_STR].append( [chrm, startPos, targetStartPos] )
                p[BREAKPOINT_STR].append([chrm, startPos, targetStartPos, DUP_STR])
                p[BREAKPOINT_STR].append([chrm, endPos, targetStartPos, DEL_STR])
            else :
                p[DUP_IMPR_STR].append( [chrm, targetStartPos, endPos] )
                p[BREAKPOINT_STR].append([chrm, targetStartPos, endPos, DUP_STR])
                p[BREAKPOINT_STR].append([chrm, targetStartPos, startPos, DEL_STR])
        elif info["SVTYPE"] == BND_STR :
            if (varID in traBuf) == False :
                traBuf[varID] = []
            traBuf[varID].append(s)

            if len(traBuf[varID]) == 6 :
                chrm = traBuf[varID][1][0]
                startPos = int(traBuf[varID][1][1])
                endPos = int(traBuf[varID][3][1])
                varLen = (endPos - startPos) + 1

                targetChrm = traBuf[varID][5][0]
                targetStartPos = int(traBuf[varID][5][1])
                targetEndPos = targetStartPos + varLen - 1
                p[TRA_STR].append( [chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )
                p[BREAKPOINT_STR].append([chrm, startPos, endPos, DEL_STR])
                if startPos < targetStartPos :
                    p[DUP_IMPR_STR].append( [chrm, startPos, targetStartPos] )
                    p[BREAKPOINT_STR].append([chrm, startPos, targetStartPos, DUP_STR])
                    p[BREAKPOINT_STR].append([chrm, endPos, targetStartPos, DEL_STR])
                else :
                    p[DUP_IMPR_STR].append( [chrm, targetStartPos, endPos] )
                    p[BREAKPOINT_STR].append([chrm, targetStartPos, endPos, DUP_STR])
                    p[BREAKPOINT_STR].append([chrm, targetStartPos, startPos, DEL_STR])
    p_file.close()

    '''
    ff = open("./tp.txt","w")
    for v in p[DEL_STR] :
        ff.write("DEL\t%d\t%d\t%d\n" % (0, v[1], v[2]))
    for v in p[INV_STR] :
        ff.write("INV\t%d\t%d\t%d\n" % (0, v[1], v[2]))
    for v in p[DUP_STR] :
        ff.write("DUP\t%d\t%d\t%d\n" % (0, v[1], v[2]))
    for v in p[TRA_STR] :
        startPos = v[1]
        endPos = v[2]
        targetStartPos = v[4]
        if startPos < targetStartPos :
            ff.write("TRA\t%d\t%d\t%d\n" % (0, startPos, targetStartPos))
        else :
            ff.write("TRA\t%d\t%d\t%d\n" % (0, targetStartPos, endPos))
    ff.close()
    '''

    return p

def dgv(param) :

    # init.
    p = {}
    p[BREAKPOINT_STR] = []
    p[DEL_STR] = []
    p[DUP_STR] = []

    # file loading
    p_name = param[0]
    p_file = open(c.CONFIG["p-dir"] + "/" + p_name, "r")

    pl_microarray = ["Affymetrix_6.0", "IlluminaOmni1-Quad", "Custom_Illumina_1.2M", "Affymetrix6.0_Illumina1M", "Affymetrix6.0", "NimbleGen42M", "CytoScanHD_2.7M", "AgilentCustom_015685+015686+244K", "Illumina1M", "Agilent24M", "AgilentCustom180K_v2.1+v3.0"]
    pl_sequencing = ["SOLiD", "Illumina_II_IIX_HiSeq", "Illumina_HiSeq", "Illumina_GA+454", "SangerSequencing", "Illumina_GA", "Multiple_NGS_Sequencing"]

    considered_platforms = []
    if len(param) > 1:
        if ("seq" in param[1]) != False :
            considered_platforms += pl_sequencing
        elif "array" in param[1] :
            considered_platforms += pl_microarray
    else :
        considered_platforms += pl_microarray
        considered_platforms += pl_sequencing

    print considered_platforms

    for line in p_file :
        if line[0] == "#" :
            continue

        s = line.strip().split()

        chrm = s[0]
        info = parseVcfInfo(s[8]) 
        variant_type = info["variant_type"]
        platforms = info["Platforms"].split(",")

        # filter by platforms
        skipThis = True
        for pl in platforms :
            if (pl in considered_platforms) :
                skipThis = False
                break
        if skipThis :
            continue

        #  deletion
        if info["variant_sub_type"] == "Loss" :
            startPos = int(info["inner_start"]) # high-confident region
            endPos = int(info["inner_end"])
            #startPos = [ int(info["outer_start"]), int(info["inner_start"]) ]
            #endPos = [ int(info["inner_end"]), int(info["outer_end"]) ]
            p[DEL_STR].append([chrm, startPos, endPos])
            p[BREAKPOINT_STR].append([chrm, startPos, endPos, DEL_STR])

        '''
        if info["variant_sub_type"] == "Gain" :
            startPos = int(info["inner_start"]) # high-confident region
            endPos = int(info["inner_end"])
            #startPos = [ int(info["outer_start"]), int(info["inner_start"]) ]
            #endPos = [ int(info["inner_end"]), int(info["outer_end"]) ]
            #p[DUP_STR].append([chrm, startPos, endPos])
            #p[BREAKPOINT_STR].append([chrm, startPos, endPos, DUP_STR])
        '''

    print len(p[DEL_STR])
    return p

def giab(param) :

    # init.
    p = {}
    p[BREAKPOINT_STR] = []
    p[DEL_STR] = []
    p[INS_STR] = []

    # file loading (deletion)
    p_name = param[0]
    p_file = open(c.CONFIG["p-dir"] + "/" + p_name, "r")
    for line in p_file :
        if "Start" in line :
            continue

        s = line.strip().split()

        chrm = "chr" + s[0]
        startPos = int(s[1])
        endPos = int(s[2])
        svSize = endPos - startPos + 1
        p[DEL_STR].append([chrm, startPos, endPos])
        p[BREAKPOINT_STR].append([chrm, startPos, endPos, DEL_STR]) 

    # file loading (de novo insertion)
    '''
    p[INS_STR] = []
    p_name = param[0].replace(".bed", "_ins.bed")
    p_file = open(c.CONFIG["p-dir"] + "/" + p_name, "r")
    for line in p_file :
        if "Start" in line :
            continue

        s = line.strip().split()

        chrm = "chr" + s[0]
        startPos = int(s[1])
        endPos = int(s[2])
        svSize = endPos - startPos + 1
        p[INS_STR].append([chrm, startPos, endPos])
        p[BREAKPOINT_STR].append([chrm, startPos, endPos, INS_STR])
    '''
    print len(p[DEL_STR])
    return p


def use_overlap(param) :
    print "Use overlap"
    print param
    #indir = c.CONFIG["result-dir"] + "/" + param[1]
    #reads, aligners, vcs, vc_time = getListFromDir(indir)

    p = {}
    p["overlap"] = True
    return p


####################################
# For each variation caller
####################################
def crest(param) :
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []
    result[ITX_STR] = []

    infile = open(param[1] + "/out.sorted.bam.predSV.txt", "r")
    for line in infile :
        if line[0] == "#" :
            continue

        #0      1           2   3   4       5           6   7  8
        #chr22  16909016    +   3   chr22   16910229    +   1  DEL
        s = line.strip().split()       
        chrm = s[0]
        leftPos = int(s[1])
        rightPos = int(s[5])
        supp = int(s[3]) + int(s[7])
        svType = s[8]

        if leftPos < rightPos :
            startPos = leftPos
            endPos = rightPos
        else :
            startPos = rightPos
            endPos = leftPos

        if svType == DEL_STR :
            result[DEL_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DEL_STR])
        if svType == INS_STR :
            # potential duplications
            result[DUP_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DUP_STR])
        elif svType == INV_STR :
            result[INV_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, INV_STR])
        elif svType == ITX_STR  :
            result[INV_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, INV_STR])
    infile.close()

    adjacency = param[3]
    result = findImprDup(result)
    result = findDuplication(result, adjacency)
    result = findTranslocation(result, adjacency)

    return result

def gasvpro(param) :
    # runtime & memory
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []

    # variations
    if "150r1." in param[1] : # not working without source mod.
        return result

    if "200-20i1" in param[1] : # not working without source mod.
        return result


    infile = open(param[1] + "/BAMToGASV_AMBIG.gasv.combined.in.clusters.GASVPro.clusters", "r")
    for line in infile :
        if line[0] == "#" :
            continue

        s = line.strip().split()
        chrm = "chr" + s[1]
        leftBegin = int(s[2].split(",")[0])
        leftEnd = int(s[2].split(",")[1])
        rightBegin = int(s[4].split(",")[0])
        rightEnd = int(s[4].split(",")[1])
        supp = int(s[5])
        svType = s[7]

        leftPos = [leftBegin, leftEnd]
        rightPos = [rightBegin, rightEnd]

        if svType == "D" :
            result[DEL_STR].append( [supp, chrm, leftPos, rightPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, leftPos, rightPos, DEL_STR])
        elif svType[0] == "I" :
            result[INV_STR].append( [supp, chrm, leftPos, rightPos] )
            result[BREAKPOINT_STR].append([supp, chrm, leftPos, rightPos, INV_STR])
    infile.close()

    return result

def lumpy(param) :
    # runtime & memory
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []

    infile = open(param[1] + "/out.vcf", "r")
    for line in infile :
        if line[0] == "#" :
            continue

        #does not provide this information
        #if s[6] != "PASS" :
        #    continue
        s = line.strip().split()

        chrm = s[0]
        pos = int(s[1])
        varID = s[2]
        info = parseVcfInfo(s[7]) 
        supp = int(float(info["SU"]))

        if info["SVTYPE"] == DEL_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[DEL_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DEL_STR])
        if info["SVTYPE"] == INS_STR :
            startPos = pos 
            endPos = pos + int(info["SVLEN"]) - 1
            result[INS_STR].append( [supp, chrm, startPos, endPos ] )
        elif info["SVTYPE"] == INV_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[INV_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, INV_STR])
        elif info["SVTYPE"] == DUP_STR  :
            startPos = pos
            endPos = int(info["END"])
            result[DUP_STR].append( [supp, chrm, startPos, endPos] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DUP_STR])
    infile.close()

    adjacency = param[3]
    result = findImprDup(result)
    result = findDuplication(result, adjacency)
    result = findTranslocation(result, adjacency)

    return result

def gustaf(param) :
    # runtime & memory
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []

    traBuf = {}
    infile = open(param[1] + "/out.vcf", "r")
    for line in infile :
        if line[0] == "#" :
            continue

        s = line.strip().split()

        chrm = s[0]
        pos = int(s[1])
        varID = s[2]
        info = parseVcfInfo(s[7]) 
        supp = int(info["DP"])

        if ("SVTYPE" in info) == False :
            info["SVTYPE"] = "SNP"
        if info["SVTYPE"] == DEL_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[DEL_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DEL_STR])
        if info["SVTYPE"] == INS_STR :
            startPos = pos 
            endPos = pos + int(info["SVLEN"]) - 1
            result[INS_STR].append( [supp, chrm, startPos, endPos ] )
        elif info["SVTYPE"] == INV_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[INV_STR].append( [supp, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, INV_STR])
        elif info["SVTYPE"] == DUP_STR :
            if ("IMPRECISE" in info) or ":TANDEM" in s[4]  :
                startPos = pos
                endPos = int(info["END"])
                result[DUP_IMPR_STR].append( [supp, chrm, startPos, endPos] )
                result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DUP_STR])
            else :
                startPos = pos 
                endPos = pos + abs(int(info["SVLEN"])) - 1
                target = int(info["TARGETPOS"])
                targetChrm = chrm
                targetStartPos = target
                targetEndPos = targetStartPos + abs(int(info["SVLEN"])) - 1
                result[DUP_STR].append( [supp, chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )

                if startPos < targetStartPos :
                    result[BREAKPOINT_STR].append([supp, chrm, startPos, targetStartPos, DUP_STR])
                    result[BREAKPOINT_STR].append([supp, chrm, endPos, targetStartPos, DEL_STR])
                else :
                    result[BREAKPOINT_STR].append([supp, chrm, targetStartPos, endPos, DUP_STR])
                    result[BREAKPOINT_STR].append([supp, chrm, targetStartPos, startPos, DEL_STR])
        elif info["SVTYPE"] == TRA_STR :
            if varID == "." : # SNPs
                continue
            _s = varID.split("_")
            varID = _s[0] + "_" + _s[1]

            if (varID in traBuf) == False :
                traBuf[varID] = []
            traBuf[varID].append( s + [supp] )

            if len(traBuf[varID]) == 6 :
                chrm = traBuf[varID][1][0]
                startPos = int(traBuf[varID][1][1])
                endPos = int(traBuf[varID][3][1])
                varLen = (endPos - startPos) + 1

                targetChrm = traBuf[varID][5][0]
                targetStartPos = int(traBuf[varID][5][1])
                targetEndPos = targetStartPos + varLen - 1

                suppList = []
                for i in range(0, 6) :
                    suppList.append( traBuf[varID][i][-1] )
                minSupp = min( suppList )
                result[TRA_STR].append( [minSupp, chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )

                result[BREAKPOINT_STR].append([minSupp, chrm, startPos, endPos, DEL_STR])
                if startPos < targetStartPos :
                    result[BREAKPOINT_STR].append([minSupp, chrm, startPos, targetStartPos, DUP_STR])
                    result[BREAKPOINT_STR].append([minSupp, chrm, endPos, targetStartPos, DEL_STR])
                else :
                    result[BREAKPOINT_STR].append([minSupp, chrm, targetStartPos, endPos, DUP_STR])
                    result[BREAKPOINT_STR].append([minSupp, chrm, targetStartPos, startPos, DEL_STR])
    infile.close()
    
    return result 

def v(param) :
    indir = param[1]

    result = {}
    for t in VARIATIONS :
        result[t] = []

    # deletion
    #infile = open(param[1] + "/del_KL.txt", "r")
    infile = open(param[1] + "/del.txt", "r")

    for line in infile :
        s = line.strip().split()
        chrm = s[0]
        pos1 = int(s[2])
        pos2 = int(s[3])
        supp = int(float(s[4]))

        result[DEL_STR].append( [supp, chrm, pos1, pos2] ) # [start, end], 1-based
    infile.close()

    infile = open(param[1] + "/time.log", "r")
    lastLine =  ""
    for line in infile :
        lastLine = line
    infile.close()
    timeSpan = int(lastLine)
    result[RUNTIME_STR] = timeSpan

    return result

def v_seqan(param) :
    #param = [p, path, overlap, adjacency, table]

    result = getRuntimeAndMemory(param[1] + "/log.txt")
    for t in VARIATIONS :
        result[t] = []

    traBuf = {}
    infile = open(param[1] + "/out.vcf", "r")
    outfile = open("./fp.txt", "w")
    for line in infile :
        if line[0] == "#" :
            continue
        s = line.strip().split()

        chrm = s[0]
        pos = int(s[1])
        varID = s[2]
        info = parseVcfInfo(s[7]) 
        rd = float(info["RD"])
        vt = float(info["VT"])
        gc = float(info["GC"])
        cp = float(info["CP"])
        se = int(info["SE"])
        pe = int(info["PE"])
        ce = int(info["CE"])
        re = float(info["RE"])

        # ce : run to get soft-clipped reads related stat. (no filtering)
        if param[4] != "ce" and s[6] != "PASS" :
            continue

        if s[6] == "MERGED" :
            continue
        
        if "RT" in info : # old ver.
            supp = int(float(info["RT"]))
        elif param[4] == "vt" :
            supp = int(vt);
        else :
            supp = float(info["SC"])
        
        suppList = [supp, se, pe, ce, re, vt, gc, cp, rd]

        if info["SVTYPE"] == DEL_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[DEL_STR].append( [suppList, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append( [suppList, chrm, startPos, endPos, DEL_STR] )

        elif info["SVTYPE"] == INS_STR :
            startPos = pos 
            endPos = pos + int(info["SVLEN"]) - 1
            result[INS_STR].append( [suppList, chrm, startPos, endPos ] )

        elif info["SVTYPE"] == INV_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            result[INV_STR].append( [suppList, chrm, startPos, endPos ] )
            result[BREAKPOINT_STR].append([suppList, chrm, startPos, endPos, INV_STR])

        elif info["SVTYPE"] == DUP_STR :
            if ("IMPRECISE" in info) or ":TANDEM" in s[4]  :
                startPos = pos
                endPos = pos + abs(int(info["SVLEN"])) - 1
                result[DUP_IMPR_STR].append( [suppList, chrm, startPos, endPos] )
                result[BREAKPOINT_STR].append([suppList, chrm, startPos, endPos, DUP_STR])
            else :
                startPos = pos 
                endPos = pos + abs(int(info["SVLEN"])) - 1
                target = int(info["TARGETPOS"])
                targetChrm = chrm
                targetStartPos = target
                targetEndPos = targetStartPos + abs(int(info["SVLEN"])) - 1
                result[DUP_STR].append( [suppList, chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )

                if startPos < targetStartPos :
                    result[DUP_IMPR_STR].append( [suppList, chrm, startPos, targetStartPos] )
                    result[BREAKPOINT_STR].append([suppList, chrm, startPos, targetStartPos, DUP_STR])
                    result[BREAKPOINT_STR].append([suppList, chrm, endPos, targetStartPos, DEL_STR])
                else :
                    result[DUP_IMPR_STR].append( [suppList, chrm, targetStartPos, endPos] )
                    result[BREAKPOINT_STR].append([suppList, chrm, targetStartPos, endPos, DUP_STR])
                    result[BREAKPOINT_STR].append([suppList, chrm, targetStartPos, startPos, DEL_STR])

        elif info["SVTYPE"] == TRA_STR :
            startPos = pos 
            endPos = pos + abs(int(info["SVLEN"])) - 1
            target = int(info["TARGETPOS"])
            targetChrm = chrm
            targetStartPos = target
            targetEndPos = targetStartPos + abs(int(info["SVLEN"])) - 1
            result[TRA_STR].append( [suppList, chrm, startPos, endPos, targetChrm, targetStartPos, targetEndPos ] )
            result[BREAKPOINT_STR].append([suppList, chrm, startPos, endPos, DEL_STR])
             
            if startPos < targetStartPos :
                result[DUP_IMPR_STR].append( [suppList, chrm, startPos, targetStartPos] )
                result[BREAKPOINT_STR].append([suppList, chrm, startPos, targetStartPos, DUP_STR])
                result[BREAKPOINT_STR].append([suppList, chrm, endPos, targetStartPos, DEL_STR])
            else :
                result[DUP_IMPR_STR].append( [suppList, chrm, targetStartPos, endPos] )
                result[BREAKPOINT_STR].append([suppList, chrm, targetStartPos, endPos, DUP_STR])
                result[BREAKPOINT_STR].append([suppList, chrm, targetStartPos, startPos, DEL_STR])
    infile.close()

    return result 

def pindel(param) :
    # runtime & memory
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []

    # deletion
    infile = open(param[1] + "/out.txt_D", "r")
    for line in infile :
        if not line[0].isdigit() :
            continue
        s = line.strip().split()
        chrm = s[7]
        pos1 = int(s[9])
        pos2 = int(s[10])
        supp = int(s[15])
        result[DEL_STR].append( [supp, chrm, pos1, pos2] ) # [start, end], 1-based
        result[BREAKPOINT_STR].append([supp, chrm, pos1, pos2, DEL_STR])
    infile.close()

    # insertion
    infile = open(param[1] + "/out.txt_SI", "r")
    for line in infile :
        if not line[0].isdigit() :
            continue
        s = line.strip().split()
        chrm = s[7]
        pos1 = int(s[9])
        pos2 = int(s[10])
        supp = int(s[15])
        result[INS_STR].append( [supp, chrm, pos1, pos2] )
    infile.close()
    infile = open(param[1] + "/out.txt_LI", "r")
    for line in infile :
        if not line[0].isdigit() :
            continue
        s = line.strip().split()
        chrm = s[7]
        pos1 = int(s[9])
        pos2 = int(s[10])
        supp = int(s[15])
        result[INS_STR].append( [supp, chrm, pos1, pos2] ) # [start, end], 1-based
    infile.close()

    # inversion
    infile = open(param[1] + "/out.txt_INV", "r")
    for line in infile :
        if not line[0].isdigit() :
            continue
        s = line.strip().split()
        chrm = s[7]
        pos1 = int(s[9])
        pos2 = int(s[10])
        supp = int(s[15])
        result[INV_STR].append( [supp, chrm, pos1, pos2] ) # [start, end], 1-based
        result[BREAKPOINT_STR].append([supp, chrm, pos1, pos2, INV_STR])
    infile.close()

    # duplication
    infile = open(param[1] + "/out.txt_TD", "r")
    for line in infile :
        if not line[0].isdigit() :
            continue
        s = line.strip().split()
        chrm = s[7]
        pos1 = int(s[9])
        pos2 = int(s[10])
        supp = int(s[15])
        result[DUP_STR].append( [supp, chrm, pos1, pos2] ) # [start, end], 1-based
        result[BREAKPOINT_STR].append([supp, chrm, pos1, pos2, DUP_STR])
    infile.close()

    infile = open(param[1] + "/log.txt", "r")
    timeSpan = 0
    for line in infile :
        if line.strip()[-8:] != "seconds.":
            continue

        sec = line.split(" seconds")[0].split(" ")[-1]
        timeSpan += int(sec)
    infile.close()
    result[TIME_ELAPSED_STR] = timeSpan
    
    adjacency = param[3]
    result = findImprDup(result)
    result = findDuplication(result, adjacency)
    result = findTranslocation(result, adjacency)

    return result

def _delly(param) :
    # variations
    result = {}
    for t in VARIATIONS :
        result[t] = []

    if os.path.isfile(param[1] + "/out.vcf") :
        infile = open(param[1] + "/out.vcf", "r")
        for line in infile :
            if line[0] == '#':
                continue
            s = line.strip().split()

            if s[6] != "PASS" :
                continue

            chrm = s[0]
            pos = int(s[1])
            info = parseVcfInfo(s[7]) 

            # evidence
            suppPE = 0
            suppSR = 0
            if "PE" in info :
                suppPE = int(info["PE"])
            if "SR" in info :
                suppSR = int(info["SR"])
            supp = suppPE + suppSR

            if (info["SVTYPE"] == "DEL") :
                startPos = pos
                endPos = int(info["END"])
                result[DEL_STR].append( [supp, chrm, startPos, endPos ] )
                result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DEL_STR])
            elif info["SVTYPE"] == INS_STR :
                startPos = pos 
                endPos = pos + int(info["INSLEN"]) - 1
                result[INS_STR].append( [supp, chrm, startPos, endPos ] )
            elif info["SVTYPE"] == INV_STR :
                startPos = pos 
                endPos = abs(int(info["END"]))
                result[INV_STR].append( [supp, chrm, startPos, endPos ] )
                result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, INV_STR])
            elif info["SVTYPE"] == DUP_STR :
                    startPos = pos
                    endPos = int(info["END"])
                    result[DUP_STR].append( [supp, chrm, startPos, endPos] )
                    result[DEL_STR].append( [supp, chrm, startPos, endPos] )
                    result[BREAKPOINT_STR].append([supp, chrm, startPos, endPos, DUP_STR])
        infile.close()

    return result

def delly(param) :
    # runtime & memory
    result = getRuntimeAndMemory(param[1] + "/log.txt")

    # variations
    for t in VARIATIONS :
        result[t] = []

    # merge all the result in "DEL" section
    if ("DEL" in param[1]) == False :
        return result

    _result = _delly(param)
    result[DEL_STR] = _result[DEL_STR]
    result[BREAKPOINT_STR] = _result[BREAKPOINT_STR]

    _param = param
    _param[1] = _param[1].replace("DEL", "INV")
    _result = _delly(param)
    result[INV_STR] = _result[INV_STR]
    result[BREAKPOINT_STR] += _result[BREAKPOINT_STR]

    _param = param
    _param[1] = _param[1].replace("INV", "INS")
    _result = _delly(param)
    result[INS_STR] = _result[INS_STR]
    result[BREAKPOINT_STR] += _result[BREAKPOINT_STR]

    _param = param
    _param[1] = _param[1].replace("INS", "DUP")
    _result = _delly(param)
    result[DUP_STR] = _result[DUP_STR]
    result[BREAKPOINT_STR] += _result[BREAKPOINT_STR]

    adjacency = param[3]
    result = findImprDup(result)
    result = findDuplication(result, adjacency)
    result = findTranslocation(result, adjacency)
    
    print len(result[TRA_STR])
    #print (result[TRA_STR])
    #exit(1)

    return result


####################################
def isAdjacent(beginPos1, endPos1, beginPos2, endPos2, tol) :
    return ( (endPos1 <= beginPos2 and ((beginPos2 - endPos1) <= tol)) or \
             (endPos2 <= beginPos1 and ((beginPos1 - endPos2) <= tol)) );

def findImprDup(result) :
    if len(result[DUP_IMPR_STR]) == 0 :
        for sv in result[DUP_STR] :
            if len(sv) == 4 :
                result[DUP_IMPR_STR].append(sv)
            else :
                # sv[0~4] : supp, chrm, begin, end, targetChrm, targetBegin
                result[DUP_IMPR_STR].append( [sv[0], sv[1], sv[2], sv[5]] )
    return result


def findDuplication4Gasv(result, diff) :
    result[DUP_STR] = []

    listCheckRemoval = []
    for i in xrange(0, len(result[DEL_STR])) :
        listCheckRemoval.append(False)

    # check DEL
    for sv_del in result[DEL_STR] :

        leftSupp = sv_del[0]
        rightSupp = sv_del[0]

        for sv_impr_dup in list(result[DUP_IMPR_STR]) :
            # sv_impr_dup[0~3] : supp, chrm, begin, targetPos
            # sv_del[0~3] : supp, chrm, begin, end

            # left overlap
            if (sv_impr_dup[1] == sv_del[1]) and (isMatched(sv_impr_dup[2], sv_del[2], diff)) :
                chrm = sv_del[1]
                #supp = sv_impr_dup[0] + sv_del[0]
                supp = sv_impr_dup[0]
                size = sv_impr_dup[3] - sv_del[3]
                result[DUP_STR].append( [supp, chrm, sv_del[3], sv_impr_dup[3], chrm, sv_impr_dup[2], sv_impr_dup[2] - size] )
                result[DUP_IMPR_STR].remove(sv_impr_dup)

            # right overlap
            elif (sv_impr_dup[1] == sv_del[1]) and (isMatched(sv_impr_dup[3], sv_del[3], diff)) :
                chrm = sv_del[1]
                supp = sv_impr_dup[0] + sv_del[0]
                size = sv_del[2] - sv_impr_dup[2]
                result[DUP_STR].append( [supp, chrm, sv_impr_dup[2], sv_del[2], chrm, sv_impr_dup[3], sv_impr_dup[3] + size] )
                result[DUP_IMPR_STR].remove(sv_impr_dup)

    return result

def findDuplication(result, diff) :
    # already found
    if (len(result[DUP_STR])) > 0 and (len(result[DUP_STR][0]) > 4) :
        return result

    result[DUP_STR] = []

    # check DEL
    for sv_impr_dup in list(result[DUP_IMPR_STR]) :
        for sv_del in list(result[DEL_STR]) :
            # sv_impr_dup[0~3] : supp, chrm, begin, targetPos
            # sv_del[0~3] : supp, chrm, begin, end

            # left overlap
            if (sv_impr_dup[1] == sv_del[1]) and (sv_del[3] < sv_impr_dup[3]) and (isMatched(sv_impr_dup[2], sv_del[2], diff)) :
                chrm = sv_del[1]
                #supp = sv_impr_dup[0] + sv_del[0]
                supp = sv_impr_dup[0]

                if type(sv_impr_dup[2]) is not type(0) :
                    size = sv_impr_dup[3][0] - sv_del[3][1]
                    result[DUP_STR].append( [supp, chrm, sv_del[3], sv_impr_dup[3], chrm, sv_impr_dup[2], [sv_impr_dup[2][0] - size, sv_impr_dup[2][1] - size]] )
                else :
                    size = sv_impr_dup[3] - sv_del[3]
                    result[DUP_STR].append( [supp, chrm, sv_del[3], sv_impr_dup[3], chrm, sv_impr_dup[2], sv_impr_dup[2] - size] )
                
                result[DEL_STR].remove(sv_del)

            # right overlap
            elif (sv_impr_dup[1] == sv_del[1]) and (sv_impr_dup[2] < sv_del[2]) and (isMatched(sv_impr_dup[3], sv_del[3], diff)) :
                chrm = sv_del[1]
                supp = sv_impr_dup[0] + sv_del[0]
                
                if type(sv_impr_dup[2]) is not type(0) :
                    size = sv_del[2][0] - sv_impr_dup[2][1]
                    result[DUP_STR].append( [supp, chrm, sv_impr_dup[2], sv_del[2], chrm, sv_impr_dup[3], [sv_impr_dup[3][0] + size, sv_impr_dup[3][1] + size]] ) 
                else :
                    size = sv_del[2] - sv_impr_dup[2]
                    result[DUP_STR].append( [supp, chrm, sv_impr_dup[2], sv_del[2], chrm, sv_impr_dup[3], sv_impr_dup[3] + size] )
                result[DEL_STR].remove(sv_del)

    return result

def isTranslocation(dup_i, dup_j, diff) :
    if dup_i[1] != dup_j[1] : # chrm
        return False

    if dup_i[2] < dup_j[2] :
        left = dup_i
        right = dup_j
    else :
        left = dup_j
        right = dup_i

    # supp[0], chrm[1], begin[2], end[3], targetChrm[4], targetBegin[5], targetEnd[6]
    return isMatched(left[2], right[5], diff) and isMatched(right[3], left[5], diff)


def findTranslocation(result, diff) :
    # find translocations from duplications
    listCheckRemoval = []
    for i in xrange(0, len(result[DUP_STR])) :
        listCheckRemoval.append(False)

    for i in xrange(0, len(result[DUP_STR])) :
        dup_i = result[DUP_STR][i]
        for j in xrange(i+1, len(result[DUP_STR])) :
            dup_j = result[DUP_STR][j]

            if listCheckRemoval[i] == True or listCheckRemoval[j] == True :
                continue

            # find matches
            if isTranslocation(dup_i, dup_j, diff) :
                listCheckRemoval[i] = True
                listCheckRemoval[j] = True              
                if dup_i[2] < dup_j[2] :
                    dup_i[0] = dup_i[0] + dup_j[0]
                    result[TRA_STR].append(dup_i)
                else :
                    dup_j[0] = dup_i[0] + dup_j[0]
                    result[TRA_STR].append(dup_j)

    listNewDup = []
    for i in xrange(0, len(result[DUP_STR])) :
        if listCheckRemoval[i] == False : # true duplication
            listNewDup.append(result[DUP_STR][i])

    result[DUP_STR] = listNewDup
    return result

def getRuntimeAndMemory(filePath) :
    userTime = 0
    systemTime = 0
    elapsedTime = 0
    memoryUsage = 0
    if os.path.isfile(filePath) :
        infile = open(filePath, "r")
        for line in infile :
            lineStr = "User time (seconds):"
            if (lineStr in line) != False :
                userTime = float(line.split(lineStr)[1].strip())

            lineStr = "System time (seconds):"
            if (lineStr in line) != False:
                systemTime = float(line.split(lineStr)[1].strip())

            lineStr = "Elapsed (wall clock) time (h:mm:ss or m:ss):"
            if (lineStr in line) != False:
                elapsedTime = line.split(lineStr)[1].strip().split(":")
                if len(elapsedTime) == 2 :
                    elapsedTime = float(elapsedTime[0]) * 60 + float(elapsedTime[1])
                else :
                    elapsedTime = float(elapsedTime[0]) * 3600 + float(elapsedTime[1]) * 60 + float(elapsedTime[2])

            lineStr = "Maximum resident set size (kbytes):"
            if (lineStr in line) != False:
                memoryUsage = float(line.split(lineStr)[1].strip())
        infile.close()
            
    result = {}
    result[TIME_ELAPSED_STR] = elapsedTime
    result[TIME_CPU_STR] =  userTime + systemTime
    result[MEMORY_STR] = memoryUsage

    return result

####################################
# validations
####################################
def filterOutResult(vc_result, validChrm, minLen, MaxLen, minSupp = 0) :

    filtered_result = {}
    for sv_type in vc_result :

        # CPU time etc.
        if type(vc_result[sv_type]) is not type([]) :
            filtered_result[sv_type] = vc_result[sv_type]
        else :
            filtered_result[sv_type] = []
            if len(vc_result[sv_type]) == 0 :
                continue

            # get quantile
            _minSupp = minSupp
            if (type(minSupp) !=  type(0)) and minSupp.isdigit() == False : # "q"
                if type(vc_result[sv_type][0][0]) is not type(0) : 
                    suppRead = [ x[0][0] for x in vc_result[sv_type] ]
                else :
                    suppRead = [ x[0] for x in vc_result[sv_type] ] # x[0] : num of supporting reads
                suppRead.sort(reverse=True)
                q1 = suppRead[int(math.ceil(len(suppRead) * 0.25)) - 1]
                q2 = suppRead[int(math.ceil(len(suppRead) * 0.50)) - 1]
                q3 = suppRead[int(math.ceil(len(suppRead) * 0.75)) - 1]
                q4 = suppRead[len(suppRead) - 1]

                if type(vc_result[sv_type][0][0]) is not type(0) : 
                    q1_elem = [x for x in vc_result[sv_type] if x[0][0] >= q1]
                    q2_elem = [x for x in vc_result[sv_type] if x[0][0] >= q2]
                    q3_elem = [x for x in vc_result[sv_type] if x[0][0] >= q3]
                    q4_elem = [x for x in vc_result[sv_type] if x[0][0] >= q4]
                else :
                    q1_elem = [x for x in vc_result[sv_type] if x[0] >= q1]
                    q2_elem = [x for x in vc_result[sv_type] if x[0] >= q2]
                    q3_elem = [x for x in vc_result[sv_type] if x[0] >= q3]
                    q4_elem = [x for x in vc_result[sv_type] if x[0] >= q4]

                if _minSupp == "q1" :
                    _minSupp = q1
                    vc_result[sv_type] = q1_elem 
                elif _minSupp == "q2" :
                    _minSupp = q2
                    vc_result[sv_type] = q2_elem 
                elif _minSupp == "q3" :
                    _minSupp = q3
                    vc_result[sv_type] = q3_elem 
                elif _minSupp == "q4" :
                    _minSupp = q4
                    vc_result[sv_type] = q4_elem 
                else :
                    print "minSupp error : %s" % minSupp
                    exit(1)
            else :
                _minSupp = int(minSupp)

            for item in vc_result[sv_type] :
                # item : [2, 'chr10', 383179, 383254]
                # this chromosome is not considered in the comparison
                if (item[1] in validChrm) == False :
                    continue

                # get size
                if type(item[2]) is not type(0) :
                    svlen = item[3][0] - item[2][1] + 1
                    result_by_range = True
                else :
                    svlen = item[3] - item[2] + 1

                # get supp
                supp = 0
                if type(item[0]) is not type(0) :
                    supp = item[0][0] 
                else :
                    supp = item[0]

                # fiileter by size
                if (svlen < minLen) or (svlen > maxLen) or (supp < _minSupp):
                    continue
                filtered_result[sv_type].append(item)

    return filtered_result

def getListFromDir(indir) :  
    reads = []
    aligners = []
    vcs = []
    vcs_with_time = []
    for read in [x for x in os.listdir(indir)]:
        path1 = indir + "/" + read
        if not os.path.isdir(path1) :
            continue

        # read filter
        if len(c.CONFIG["read-white-list"]) > 0 :
            found = False;
            for keyword in c.CONFIG["read-white-list"] :
                if keyword in read :
                    found = True
                    break
            if found == False :
                continue

        print "Load dataset(read) : %s" % read
        if (read in reads) == False :
            reads.append(read)

        for aligner in os.listdir(path1) :
            path2 = path1 + "/" + aligner
            if not os.path.isdir(path2) :
                continue

            # mapper filter
            if len(c.CONFIG["mapper-white-list"]) > 0 :
                found = False;
                for keyword in c.CONFIG["mapper-white-list"] :
                    if keyword in aligner :
                        found = True
                        break
                if found == False :
                    continue

            if (aligner in aligners) == False :
                aligners.append(aligner)

            for vc in os.listdir(path2) :
                path3 = path2 + "/" + vc
                if not os.path.isdir(path3) :
                    continue

                if (("delly" in vc) == True) and ( (("DEL" in vc) == False) or (("q:20" in vc) == False) ) :
                    continue

                # vc filter
                if len(c.CONFIG["vc-white-list"]) > 0 :
                    found = False;
                    for keyword in c.CONFIG["vc-white-list"] :
                        if keyword in vc :
                            found = True
                            break
                    if found == False :
                        continue

                modified_time = os.stat(path3+"/log.txt").st_mtime
                vc_info = [modified_time, vc]
                vc_name = vc.split("#")[0]
                if vc_name in c.CONFIG["variant-callers"] :
                    print "Load dataset(variant caller) : %s" % vc
                    vcs_with_time.append(vc_info)

    reads.sort()
    aligners.sort()
    vcs_with_time.sort(key=lambda x:x[0], reverse=True)

    vc_time = {}
    for i in range(0, len(vcs_with_time)) :
        vc_info = vcs_with_time[i];
        vc_info_name = vc_info[1]
        vc_info_time = time.ctime(vc_info[0]);

        if (vc_info_name in vcs) == False :
            vcs.append(vc_info_name)
            vc_time[vc_info_name] = vc_info_time


    return reads, aligners, vcs, vc_time

def getMendelianError(vc_result_by_key, overlap) :

    # get union
    p = {}
    keys = vc_result_by_key.keys()
    keys.sort()
    print keys

    offstpring_key = keys[0]
    for sv_type in vc_result_by_key[offstpring_key]:
        if (sv_type in VARIATIONS) == False :
            continue
        if (sv_type in p) == False :
            p[sv_type] = []
   
        # sv[0] : evidences
        # sv[1] : chrm
        # sv[2] : beginpos
        # sv[3] : endpos
        # sv[4] : breakpoint type (BP)
        for sv in vc_result_by_key[offstpring_key][sv_type] :
            if sv_type == BREAKPOINT_STR :
                # 0 : chrm
                # 1 : sub-bp
                # 2 : begin
                # 3 : end
                # 4 : [key]
                p[sv_type].append( [sv[1], sv[4], sv[2], sv[3], [offstpring_key]] )
            else :
                p[sv_type].append( [sv[1], "", sv[2], sv[3], [offstpring_key]] )

    # find overlaps
    for key in keys[1:]:
        print key
        for sv_type in vc_result_by_key[key]:
            if (sv_type in VARIATIONS) == False :
                continue

            for sv in vc_result_by_key[key][sv_type] :
                found = False

                # sv : [[76, 2, 4, 0, 0.732714, 3.0, 0.423, 1.981084, 3.414286], 'chr22', 40721213, 40725473]
                # sv[0] : evidences
                # sv[1] : chrm
                # sv[2] : beginpos
                # sv[3] : endpos
                # sv[4] : sub-bp (BP)
                #
                # target_sv[0] : chrm
                # target_sv[1] : sub-bp
                # target_sv[2] : begin
                # target_sv[3] : end
                # target_sv[4] : key
                for target_sv in range(0, len(p[sv_type])):
                    # chrm
                    if (sv[1] != p[sv_type][target_sv][0]) :
                        continue

                    # sub bp type
                    if (sv_type == BREAKPOINT_STR) and sv[4] != p[sv_type][target_sv][1] :
                        continue

                    targetBeginPos = p[sv_type][target_sv][2]
                    targetEndPos = p[sv_type][target_sv][3]

                    if isReciOverlap(sv[2], sv[3], targetBeginPos, targetEndPos, overlap) == True :
                    #if isOverlap(sv[2], sv[3], targetBeginPos, targetEndPos) == True :
                        p[sv_type][target_sv][-1].append(key)
                        #break
    return p

def getUnion(vc_result_by_key, overlap) :

    # get union
    p = {}
    keys = vc_result_by_key.keys()
    keys.sort()
    print keys

    for key in keys:
        for sv_type in vc_result_by_key[key]:
            if (sv_type in VARIATIONS) == False :
                continue
            if (sv_type in p) == False :
                p[sv_type] = []
   
            # find union
            for new_sv in vc_result_by_key[key][sv_type] :
                found = False

                # new_sv : [[76, 2, 4, 0, 0.732714, 3.0, 0.423, 1.981084, 3.414286], 'chr22', 40721213, 40725473]
                # new_sv[0] : evidences
                # new_sv[1] : chrm
                # new_sv[2] : beginpos
                # new_sv[3] : endpos
                # new_sv[4] : breakpoint type (BP)
                # find overlaps
                listOverlap = []
                for old_i in range(0, len(p[sv_type])):
                    # chrm
                    if (new_sv[1] != p[sv_type][old_i][0]) :
                        continue

                    # sub bp type
                    if (sv_type == BREAKPOINT_STR) and new_sv[4] != p[sv_type][old_i][1] :
                        continue

                    # check intervals
                    # [0] chrm
                    # [1] sub_bp_type (BREAKPOINT)
                    # [2] list of invertvals
                    # [3] supports
                    for interval in p[sv_type][old_i][2] :
                        if isReciOverlap(new_sv[2], new_sv[3], interval[0], interval[1], overlap) == True :
                            listOverlap.append(old_i)
                            break

                # new sv
                if len(listOverlap) == 0 :
                        # [0] chrm
                        # [1] sub_bp_type (BREAKPOINT)
                        # [2] list of invertvals
                        # [3] supports
                        if sv_type == BREAKPOINT_STR :
                            p[sv_type].append( [new_sv[1], new_sv[4], [[new_sv[2], new_sv[3]]], []] )
                        else :
                            p[sv_type].append( [new_sv[1], "", [[new_sv[2], new_sv[3]]], []] )
                else :
                    # copy one
                    mergedSV = copy.deepcopy(p[sv_type][listOverlap[0]])

                    # extract intervals
                    intervalList = [ [new_sv[2], new_sv[3]] ]
                    for i in listOverlap :
                        intervalList += p[sv_type][i][2]
                    mergedSV[2] = copy.deepcopy(intervalList)

                    # remove
                    for i in sorted(listOverlap, reverse=True):
                        del p[sv_type][i]

                    p[sv_type].append(mergedSV)

    # get overlaps
    for key in keys:
        for sv_type in vc_result_by_key[key]:
            if (sv_type in VARIATIONS) == False :
                continue

            # find overlaps
            for sv in vc_result_by_key[key][sv_type] :
                found = False

                # new_sv : [[76, 2, 4, 0, 0.732714, 3.0, 0.423, 1.981084, 3.414286], 'chr22', 40721213, 40725473]
                # new_sv[0] : evidences
                # new_sv[1] : chrm
                # new_sv[2] : beginpos
                # new_sv[3] : endpos
                # new_sv[4] : breakpoint type (BP)
                for union_sv in range(0, len(p[sv_type])):
                    # chrm
                    if (sv[1] != p[sv_type][union_sv][0]) :
                        continue

                    # sub bp type
                    if (sv_type == BREAKPOINT_STR) and sv[4] != p[sv_type][union_sv][1] :
                        continue

                    for interval in p[sv_type][union_sv][2] :
                        if isReciOverlap(sv[2], sv[3], interval[0], interval[1], overlap) == True :
                            p[sv_type][union_sv][-1].append(key)
                            break

    return p

def validation(p, genomeDir, overlap, adjacency, minLen, maxLen, maxTH, reportName, table, separation, minSupp) :
    result = {}
    indir = c.CONFIG["result-dir"] + "/" + genomeDir
    reads, aligners, vcs, vc_time = getListFromDir(indir)

    overlapReport = False
    if ("overlap" in p) == True and p["overlap"] == True :
        overlapReport = True

    # parsing vc result
    vc_result_by_key = {}
    print "===== Start parsing =====" 
    for read in reads :
        print read
        result[read] = {}
        for aligner in aligners :
            print "\t" + aligner
            result[read][aligner] = {}
            for vc in vcs :
                path =  indir + "/" + read + "/" +  aligner + "/" + vc
                if not os.path.isdir(path) :
                    continue
                print "\t\t" + vc

                vc_name = vc.split("#")[0]
                param = [p, path, overlap, adjacency, table]
                vc_result = eval(vc_name)(param)
                key = "%s$$%s$$%s" % (read, aligner, vc)

                if overlapReport :
                    if "v_seqan" in key :
                        key = "0." + key
                    elif "lumpy" in key :
                        key = "1." + key
                    elif "delly" in key :
                        key = "2." + key
                

                vc_result_by_key[key] = vc_result

    # make a report
    reportDir = c.CONFIG["report-dir"] + "/" + reportName
    checkAndMake(reportDir)
    outfile = open(reportDir + "/report.txt", "w")
    print "report: %s" % (reportDir + "/report.txt")

    # overlap report
    if overlapReport :
        
        # get valid chromosomes
        validChrm = {}
        '''
        for key in range(1, 23) :
            validChrm["chr%d" % key] = True
        validChrm["chrX"] = True
        validChrm["chrY"] = True
        '''
        for key in vc_result_by_key:
            for sv_type in vc_result_by_key[key]:
                if (sv_type in VARIATIONS) == False :
                    continue
                for item in vc_result_by_key[key][sv_type] :
                    validChrm[item[1]] = True

        # filtering
        keys = vc_result_by_key.keys()
        keys.sort()
        for key in keys :
            outfile.write("%s\n" % key)

            _minSupp = minSupp
            if (("q" in minSupp) == False) and ("v_seqan" in key) == True :
                _minSupp = 0 # already filtered
            print _minSupp
            vc_result_by_key[key] = filterOutResult(vc_result_by_key[key], validChrm, minLen, maxLen, _minSupp)

        # mendelian erros
        p = getMendelianError(vc_result_by_key, overlap)

        # get union / intersections
        #p = getUnion(vc_result_by_key, overlap)

        # get valid chromosomes
        validChrm = {}
        for sv_type in p :
            for item in p[sv_type] :
                validChrm[item[0]] = True

        for sv_type in p :
            outfile.write("%s\t%d\n" % (sv_type, len(p[sv_type])))
            class_result = {}
            for sv in p[sv_type] :
                class_key = ""
                for key in keys :
                    if key in sv[-1] :
                        class_key += "T\t"
                    else :
                        class_key += "F\t"
                
                if (class_key in class_result) == False :
                    class_result[class_key] = 0
                class_result[class_key] += 1

            totalCnt = len(p[sv_type])
            for class_key in class_result :
                outfile.write("%s\t%d\t%f\n" % (class_key, class_result[class_key], (float(class_result[class_key]) / float(totalCnt))))
        outfile.close()
        return 
        
    # get result
    print "===== Get result =====" 
    validChrm = {} # get valid chromosomes
    for sv_type in p :
        for item in p[sv_type] :
            validChrm[item[0]] = True
    for read in reads :
        print read
        result[read] = {}
        for aligner in aligners :
            print "\t" + aligner
            result[read][aligner] = {}
            for vc in vcs :
                path =  indir + "/" + read + "/" +  aligner + "/" + vc
                if not os.path.isdir(path) :
                    continue
                print "\t\t" + vc
                
                key = "%s$$%s$$%s" % (read, aligner, vc)
                vc_result = vc_result_by_key[key]

                _minSupp = minSupp
                # filter out result
                if ("v_seqan" in key) == True :
                    _minSupp = 0 # already filtered
                vc_result = filterOutResult(vc_result, validChrm, minLen, maxLen, _minSupp)

                # get accuracy
                if table == "ce" :
                    result[read][aligner][vc] = getAccuracyForCE(p, vc_result, overlap, minLen, maxLen, separation)
                else :
                    result[read][aligner][vc] = getAccuracy(p, vc_result, overlap, minLen, maxLen, separation)
                result[read][aligner][vc][TIME_ELAPSED_STR] = vc_result[TIME_ELAPSED_STR]
                result[read][aligner][vc][TIME_CPU_STR] = vc_result[TIME_CPU_STR]
                result[read][aligner][vc][MEMORY_STR] = vc_result[MEMORY_STR]

    outfile.write("genome\t%s\n" % genomeDir)
    for read in reads :
        outfile.write("read\t%s\n" % read)
        for aligner in aligners:
            outfile.write("aligner\t%s\n" % aligner)

            #for t in p :
            #for t in [RUNTIME_STR, DEL_STR, INV_STR] :
            #for t in [RUNTIME_STR, INV_STR] :
            #for t in [RUNTIME_STR, INV_STR, DEL_STR] :
            #for t in [RUNTIME_STR, DEL_STR, INV_STR] :
            #for t in [RUNTIME_STR, DEL_STR] :
            #for t in [RUNTIME_STR, TRA_STR, DUP_STR] :
            #for t in [RUNTIME_STR, TRA_STR, DUP_STR, INV_STR, DEL_STR] :
            #for t in [TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR, DEL_STR, INV_STR, DUP_STR, DUP_IMPR_STR, TRA_STR] :
            for t in [TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR, BREAKPOINT_STR, DEL_STR, INV_STR, DUP_STR, TRA_STR, DUP_IMPR_STR] :
            #for t in [DUP_STR] :
                outfile.write("==== type\t%s ====\n" % t)

                # first line
                maxSuppRead = -1
                for vc in vcs :
                    if (vc in result[read][aligner]) == False :
                        continue

                    if (t in result[read][aligner][vc]) == False :
                        continue

                    if t in {TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR} :
                        resultLen = 1
                    else :
                        resultLen = len(result[read][aligner][vc][t])
                    if resultLen > maxSuppRead :
                        maxSuppRead = resultLen
                # vc 
                for vc in vcs :
                    if (vc in result[read][aligner]) == False :
                        continue

                    if (t in result[read][aligner][vc]) == False :
                        continue

                    if t in {TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR} :
                        outfile.write("[%s]%s\t%d" % (vc_time[vc],vc,result[read][aligner][vc][t]))
                    else :
                        outfile.write( "[%s] (%d) %s\n" % (vc_time[vc],len(p[t]),vc) )
                        maxLen = len(result[read][aligner][vc][t])
                        maxVal = maxLen

                        if (maxVal > maxTH + 1) :
                          maxVal = maxTH + 1

                        #for suppRead in xrange(1, maxVal) :
                        for suppRead in xrange(0, maxVal) :
                            if (t in result[read][aligner][vc]) == False :
                                continue
                            val = result[read][aligner][vc][t][suppRead]
                            #outfile.write("%d\t%s\t%s\t%s\t%d\t%d\t%d\n" % (suppRead, str(val[0])[:4], str(val[1])[:4], str(val[2])[:4], val[3], val[4], val[5]))
                            outfile.write("%.2f\t%s\t%s\t%s\t%d\t%s\n" % (val[-1], str(val[0])[:4], str(val[1])[:4], str(val[2])[:4], val[5], str(val[6])[:8]))
                    outfile.write("\n")
                outfile.write("\n\n")

    # precison-recall plot
    for read in reads :
        for aligner in aligners:
            for t in [DEL_STR] :
                for vc in result[read][aligner] :
                    filename = reportDir + "/" + read + "/"+ aligner + "/" + vc 
                    checkAndMake(filename)
                    filename += "/" + t + ".txt"
                    outfile = open(filename, "w")
                    maxLen = len(result[read][aligner][vc][t])

                    for suppRead in xrange(1, maxLen) :
                        val = result[read][aligner][vc][t][suppRead]
                        outfile.write("%d\t%.2f\t%.2f\t%.2f\n" % (suppRead, val[0], val[1], val[2]))  
                    outfile.close()

def validation_old(p, genomeDir, overlap, adjacency, minLen, maxLen, maxTH, reportName, table, separation) :
    result = {}
    indir = c.CONFIG["result-dir"] + "/" + genomeDir
    reads, aligners, vcs, vc_time = getListFromDir(indir)

    # get valid chromosomes
    validChrm = {}
    for sv_type in p :
        for item in p[sv_type] :
            validChrm[item[0]] = True

    print "Start parsing" 
    for read in reads :
        print read
        result[read] = {}
        for aligner in aligners :
            print "\t" + aligner
            result[read][aligner] = {}
            for vc in vcs :
                path =  indir + "/" + read + "/" +  aligner + "/" + vc
                if not os.path.isdir(path) :
                    continue
                print "\t\t" + vc

                vc_name = vc.split("#")[0]
                param = [p, path, overlap, adjacency, table]
                vc_result = eval(vc_name)(param)

                # filter out result
                vc_result = filterOutResult(vc_result, validChrm, minLen, maxLen)

                # get accuracy
                if table == "ce" :
                    result[read][aligner][vc] = getAccuracyForCE(p, vc_result, overlap, minLen, maxLen, separation)
                else :
                    result[read][aligner][vc] = getAccuracy(p, vc_result, overlap, minLen, maxLen, separation)
                result[read][aligner][vc][TIME_ELAPSED_STR] = vc_result[TIME_ELAPSED_STR]
                result[read][aligner][vc][TIME_CPU_STR] = vc_result[TIME_CPU_STR]
                result[read][aligner][vc][MEMORY_STR] = vc_result[MEMORY_STR]

    # make a report
    reportDir = c.CONFIG["report-dir"] + "/" + reportName
    checkAndMake(reportDir)

    outfile = open(reportDir + "/report.txt", "w")
    outfile.write("genome\t%s\n" % genomeDir)
    for read in reads :
        outfile.write("read\t%s\n" % read)
        for aligner in aligners:
            outfile.write("aligner\t%s\n" % aligner)

            #for t in p :
            #for t in [RUNTIME_STR, DEL_STR, INV_STR] :
            #for t in [RUNTIME_STR, INV_STR] :
            #for t in [RUNTIME_STR, INV_STR, DEL_STR] :
            #for t in [RUNTIME_STR, DEL_STR, INV_STR] :
            #for t in [RUNTIME_STR, DEL_STR] :
            #for t in [RUNTIME_STR, TRA_STR, DUP_STR] :
            #for t in [RUNTIME_STR, TRA_STR, DUP_STR, INV_STR, DEL_STR] :
            #for t in [TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR, DEL_STR, INV_STR, DUP_STR, DUP_IMPR_STR, TRA_STR] :
            for t in [TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR, BREAKPOINT_STR, DEL_STR, INV_STR, DUP_STR, TRA_STR, DUP_IMPR_STR] :
            #for t in [DUP_STR] :
                outfile.write("==== type\t%s ====\n" % t)

                # first line
                maxSuppRead = -1
                for vc in vcs :
                    if (vc in result[read][aligner]) == False :
                        continue

                    if (t in result[read][aligner][vc]) == False :
                        continue

                    if t in {TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR} :
                        resultLen = 1
                    else :
                        resultLen = len(result[read][aligner][vc][t])
                    if resultLen > maxSuppRead :
                        maxSuppRead = resultLen
                # vc 
                for vc in vcs :
                    if (vc in result[read][aligner]) == False :
                        continue

                    if (t in result[read][aligner][vc]) == False :
                        continue

                    if t in {TIME_ELAPSED_STR, TIME_CPU_STR, MEMORY_STR} :
                        outfile.write("[%s]%s\t%d" % (vc_time[vc],vc,result[read][aligner][vc][t]))
                    else :
                        outfile.write( "[%s] (%d) %s\n" % (vc_time[vc],len(p[t]),vc) )
                        maxLen = len(result[read][aligner][vc][t])
                        maxVal = maxLen

                        if (maxVal > maxTH + 1) :
                          maxVal = maxTH + 1

                        #for suppRead in xrange(1, maxVal) :
                        for suppRead in xrange(0, maxVal) :
                            if (t in result[read][aligner][vc]) == False :
                                continue
                            val = result[read][aligner][vc][t][suppRead]
                            #outfile.write("%d\t%s\t%s\t%s\t%d\t%d\t%d\n" % (suppRead, str(val[0])[:4], str(val[1])[:4], str(val[2])[:4], val[3], val[4], val[5]))
                            outfile.write("%.2f\t%s\t%s\t%s\t%d\t%s\n" % (val[-1], str(val[0])[:4], str(val[1])[:4], str(val[2])[:4], val[5], str(val[6])[:8]))
                    outfile.write("\n")
                outfile.write("\n\n")

    # precison-recall plot
    for read in reads :
        for aligner in aligners:
            for t in [DEL_STR] :
                for vc in result[read][aligner] :
                    filename = reportDir + "/" + read + "/"+ aligner + "/" + vc 
                    checkAndMake(filename)
                    filename += "/" + t + ".txt"
                    outfile = open(filename, "w")
                    maxLen = len(result[read][aligner][vc][t])

                    for suppRead in xrange(1, maxLen) :
                        val = result[read][aligner][vc][t][suppRead]
                        outfile.write("%d\t%.2f\t%.2f\t%.2f\n" % (suppRead, val[0], val[1], val[2]))  
                    outfile.close()



########################################################################
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
    checkAndMake(c.CONFIG["report-dir"])

def run(func, param1) :
    return eval(func)(param1)

def parseVcfInfo(infostr) :
    info = {}
    for item in infostr.split(";") :
        s = item.split("=")
        if len(s) == 2 :
            info[s[0]] = s[1]
        elif item != "." :
            info[item] = ""
    return info

def isMatched(t, p, result, overlap) :
    if t == DEL_STR :
        # 0:chromosome, 1:beginPos, 2:endPos
        pLen = p[2] - p[1] + 1
        resultLen = result[2] - result[1] + 1
        matchTol = float(pLen) * (1.0 - overlap)

        if( (p[0] == result[0]) and (abs(p[1] - result[1]) <= matchTol) and (abs(p[2] - result[2])) <= matchTol ) :
            return True
        else :
            return False

def isOverlap(beginPos1, endPos1, beginPos2, endPos2) :
    return ((beginPos1 <= beginPos2 and beginPos2 < endPos1) or \
           (beginPos2 <= beginPos1 and beginPos1 < endPos2) or \
           (beginPos1 <= beginPos2 and endPos2 < endPos1) or \
           (beginPos2 <= beginPos1 and endPos1 < endPos2))

def isReciOverlap(beginPos1, endPos1, beginPos2, endPos2, overlap) :

    # pos1 : result by variation callers
    # pos2 : correct answers

    # range vs. range
    if (type(beginPos1) is not type(0)) and (type(beginPos2) is not type(0)) :
        #print "type1"
        return isOverlap(beginPos1[0], beginPos1[1], beginPos2[0], beginPos2[1]) and \
               isOverlap(endPos1[0], endPos1[1], endPos2[0], endPos2[1])
    else : 
        # range vs. exact point (eg. gasvpro)
        if (type(beginPos1) is not type(0)) and (type(beginPos2) is type(0)) :
            #print "type2"
            return isOverlap(beginPos1[0], beginPos1[1], beginPos2, beginPos2) and \
                   isOverlap(endPos1[0], endPos1[1], endPos2, endPos2)

        elif (type(beginPos1) is type(0)) and (type(beginPos2) is not type(0)) :
            #print "type3"
            cnt = 0
            if isOverlap(beginPos1, beginPos1, beginPos2[0], beginPos2[1]):
                cnt += 1
            if isOverlap(endPos1, endPos1, endPos2[0], endPos2[1]):
                cnt += 1

            if cnt == 1 :
                print beginPos1, endPos1, beginPos2, endPos2

            return isOverlap(beginPos1, beginPos1, beginPos2[0], beginPos2[1]) and \
                       isOverlap(endPos1, endPos1, endPos2[0], endPos2[1])

        else : # exact point vs. exact point
            #print "type4"
            if isOverlap(beginPos1, endPos1, beginPos2, endPos2) :
                tpRegionSize = endPos2 - beginPos2 + 1
                resultRegionSize = endPos1 - beginPos1 + 1
                positions = sorted([ beginPos1, endPos1, beginPos2, endPos2 ])
                matchedRegionSize = positions[2] - positions[1] + 1

                # reciprocal overlaps
                return (matchedRegionSize >= (tpRegionSize * overlap)) and (matchedRegionSize >= (resultRegionSize * overlap))

def isMatched(pos1, pos2, diff) :
    if (type(pos1) is not type(0)) or (type(pos2) is not type(0)) : 
        if type(pos1) is type(0) :
            pos1 = [pos1, pos1]
        if pos2 is type(0) :
            pos2 = [pos2, pos2]
        return isOverlap(pos1[0], pos1[1], pos2[0], pos2[1])
    else :
        return abs(pos1-pos2) < diff

def isContained(pos, posList, diff) :
    for p in posList :
        if isMatched(pos, p, diff) :
            return True
        if (pos + diff) < p :
            break
    return False

def getAccuracyBySuppReads(sv_type, p, vc_result, overlap, readCnt) :
    # List of SVs called by callers
    # result_item = [chrm, begin, end, (SV_TYPE), supp]
    # result_item = [chrm, [start, end], [start, end], (SV_TYPE), supp] # gasvpro
    # supp = [supp, se, pe, ce, re, vt, gc, cp, rd] # v_seqan
    filteredResult = [] 
    for result_item in vc_result :
        filteredResult.append(result_item)

    # correct answers
    # tp_item = [chrm, start, end, (SV_TYPE)]
    tpCheck = []
    for tp_item in p :
        # tp_item[-1] : # of matches, compare to filtered_result
        # tp_item[-2] : list of calls from filtered_result 
        tpCheck.append( tp_item + [[]] + [0] )

    # comparision
    redundant = 0
    falsePositives = []
    truePositives = []
    positiveSupp = 0
    for result_item in filteredResult : # for all positive calls
        found = False 
        foundIdx = -1
        for (idx,tp_item) in enumerate(tpCheck) :
            # check reciprocal overlap
            if (result_item[0] == tp_item[0]) and isReciOverlap(result_item[1], result_item[2], tp_item[1], tp_item[2], overlap) == True :
                found = True
            if found == True and (sv_type == BREAKPOINT_STR and (result_item[3] != tp_item[3])) :
                found = False

            if found == True :
                foundIdx = idx
                break

        # check types of breakpoints (DEL, INV..)
        # result :  ['chr22', 17537928, 17542631, 'DEL', 47]
        # tp : ['chr22', 17537072, 17542632, 'DUP', [['chr22', 17537074, 17542631, 'DUP', 47]], 1]

        # true or false positive
        if found == True :
            positiveSupp += result_item[-1] # add reads that support positive call

            if tpCheck[foundIdx][-1] == 0 : # set true if it's found already
                tpCheck[foundIdx][-1] = 1
                tpCheck[foundIdx][-2].append(result_item)
                truePositives.append( result_item + [tp_item[1]] + [tp_item[2]] )
            else :
                # redundant case can be exist due to overlaps
                tpCheck[foundIdx][-1] += 1
                tpCheck[foundIdx][-2].append(result_item)
                redundant += 1 
        else :
            falsePositives.append(result_item)

    if len(truePositives) == 0 :
        return [0, 0, -1.0, 0, 0, 0, 0]
    else :
        tp = len(truePositives)
        recall = float(tp) / len(tpCheck)
        #precision =  float(tp) / (len(filteredResult) - redundant)
        precision =  float(tp) / (len(filteredResult))
        f1 = ((recall * precision)*2) / (recall + precision) 

        '''
        if readCnt == 0 : # most sensitive case
            #print sv_type, tp, redundant
        
            print sv_type
            print "false positives"
            for item in falsePositives:
                print item

            print "true positives"
            for item in truePositives:
                print item 

            print "lost"
            for item in tpCheck:
                if item[-1] != True: # failed to find
                    print item

            print "total positive calls"
            for result_item in filteredResult :
                print result_item
        '''

        # result_item = [chrm, begin, end, (SV_TYPE), supp]
        # supp = [supp, se, pe, ce, re, vt, gc, cp, rd]
        '''
        if sv_type == BREAKPOINT_STR :
            outfile = open("./bp.txt","w")
            outfile.write("type\tsv\tchrm\tbegin\tend\tse\tpe\tce\tre\tvt\tgc\tcp\trd\n")

            for item in truePositives:
                outfile.write("TP\t%s\t%s\t%d\t%d\t" % (sv_type, item[0], item[1], item[2]))
                outfile.write("%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\n" % (item[4][1],item[4][2],item[4][3],item[4][4],item[4][5],item[4][6],item[4][7],item[4][8]))

            for item in falsePositives:
                outfile.write("FP\t%s\t%s\t%d\t%d\t" % (sv_type, item[0], item[1], item[2]))
                outfile.write("%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\n" % (item[4][1],item[4][2],item[4][3],item[4][4],item[4][5],item[4][6],item[4][7],item[4][8]))

            outfile.close()
        '''
        
        #return [recall, precision, f1, tp, len(tpCheck), (len(filteredResult) - redundant)]
        #return [recall, precision, f1, tp, len(tpCheck), (len(filteredResult))]
        #return [recall, precision, f1, tp, len(tpCheck), (len(filteredResult) - redundant), float(del_bp) / float((len(filteredResult) - redundant))]
        return [recall, precision, f1, tp, len(tpCheck), (len(filteredResult) - redundant), positiveSupp]

def getAccuracyByScore(p, vc_result, overlap, minLen, maxLen) :
    acc = {}
    for sv_type in p : # for each SV type
        if ((sv_type in vc_result) == True) and len(vc_result[sv_type]) > 0 :
            
            # v_seqan contains [supp, se, pe, ce, re]
            for i in range(0, len(vc_result[sv_type])) :
                if type(vc_result[sv_type][i][0]) is not type(0) :
                    readCnt = vc_result[sv_type][i][0][0]
                    vc_result[sv_type][i].append(vc_result[sv_type][i][0])
                    vc_result[sv_type][i][0] = readCnt

            acc[sv_type] = []
            suppRead = [ x[0] for x in vc_result[sv_type] ] # x[0] : num of supporting reads
            suppReadMin = int(min(suppRead))
            suppReadMax = int(min( max(suppRead), MAX_SUPP ))
            #acc[sv_type].append( [-1, -1, -1] ) # supporting read : 0
            #for readCnt in xrange(1, suppReadMax+1) :
            for readCnt in xrange(0, suppReadMax+1) :
                #resultBySuppRead = [x[1:] for x in result[t] if x[0] >= readCnt]
                resultBySuppRead = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= readCnt]
                acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], resultBySuppRead, overlap, readCnt) + [readCnt] ) 
        else : 
            acc[sv_type] = {} 
    return acc 

def getAccuracyByQuantile(p, vc_result, overlap, minLen, maxLen) :
    acc = {}
    for sv_type in p : # for each SV type
        if ((sv_type in vc_result) == True) and len(vc_result[sv_type]) > 0 :
            # v_seqan contains [supp, se, pe, ce, re]
            for i in range(0, len(vc_result[sv_type])) :
                if type(vc_result[sv_type][i][0]) is not type(0) :
                    vc_result[sv_type][i][0] = vc_result[sv_type][i][0][0]

            acc[sv_type] = []
            suppRead = [ x[0] for x in vc_result[sv_type] ] # x[0] : num of supporting reads
            suppReadMin = min(suppRead)
            suppReadMax = min( max(suppRead), MAX_SUPP )

            # get quantile
            suppRead.sort(reverse=True)
            q1 = suppRead[int(math.ceil(len(suppRead) * 0.25)) - 1]
            q2 = suppRead[int(math.ceil(len(suppRead) * 0.50)) - 1]
            q3 = suppRead[int(math.ceil(len(suppRead) * 0.75)) - 1]
            q4 = suppRead[len(suppRead) - 1]

            q1_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= q1]
            q2_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if (x[0] < q1 and x[0] >= q2)]
            q3_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if (x[0] < q2 and x[0] >= q3)]
            q4_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if (x[0] < q3 and x[0] >= q4)]


            elem = [x[1:]+[x[0]] for x in vc_result[sv_type]]
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], elem, overlap, 0)  + [0])

            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q1_elem, overlap, q1) + [q1] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q2_elem, overlap, q2) + [q2] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q3_elem, overlap, q3) + [q3] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q4_elem, overlap, q4) + [q4] )

        else : 
            acc[sv_type] = {} 
    return acc


def getAccuracyByQuantileCum(p, vc_result, overlap, minLen, maxLen) :
    acc = {}
    for sv_type in p : # for each SV type
        if ((sv_type in vc_result) == True) and len(vc_result[sv_type]) > 0 :
            # v_seqan contains [supp, se, pe, ce, re]
            for i in range(0, len(vc_result[sv_type])) :
                if type(vc_result[sv_type][i][0]) is not type(0) :
                    vc_result[sv_type][i][0] = vc_result[sv_type][i][0][0]

            acc[sv_type] = []
            suppRead = [ x[0] for x in vc_result[sv_type] ] # x[0] : num of supporting reads
            suppReadMin = min(suppRead)
            suppReadMax = min( max(suppRead), MAX_SUPP )

            # get quantile
            suppRead.sort(reverse=True)
            q1 = suppRead[int(math.ceil(len(suppRead) * 0.25)) - 1]
            q2 = suppRead[int(math.ceil(len(suppRead) * 0.50)) - 1]
            q3 = suppRead[int(math.ceil(len(suppRead) * 0.75)) - 1]
            q4 = suppRead[len(suppRead) - 1]

            q1_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= q1]
            q2_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= q2]
            q3_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= q3]
            q4_elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if x[0] >= q4]


            elem = [x[1:]+[x[0]] for x in vc_result[sv_type]]
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], elem, overlap, 0)  + [0])

            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q1_elem, overlap, q1) + [q1] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q2_elem, overlap, q2) + [q2] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q3_elem, overlap, q3) + [q3] )
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], q4_elem, overlap, q4) + [q4] )

        else : 
            acc[sv_type] = {} 
    return acc

def getAccuracyByBin(p, vc_result, overlap, minLen, maxLen) :
    acc = {}
    for sv_type in p : # for each SV type
        if ((sv_type in vc_result) == True) and len(vc_result[sv_type]) > 0 :
            # v_seqan contains [supp, se, pe, ce, re]
            for i in range(0, len(vc_result[sv_type])) :
                if type(vc_result[sv_type][i][0]) is not type(0) :
                    vc_result[sv_type][i][0] = vc_result[sv_type][i][0][0]

            acc[sv_type] = []
            suppRead = [ x[0] for x in vc_result[sv_type] ] # x[0] : num of supporting reads
            suppReadMin = min(suppRead)
            suppReadMax = min( max(suppRead), MAX_SUPP )

            elem = [x[1:]+[x[0]] for x in vc_result[sv_type]]
            acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], elem, overlap, 0) + [0] )

            # get quantile
            suppRead.sort(reverse=True)
            for i in range(1, 11):
                start_th_idx = int(math.ceil(len(suppRead) * (0.1*(i-1)))) - 1
                end_th_idx = int(math.ceil(len(suppRead) * (0.1*i))) - 1

                if start_th_idx < 0 :
                    max_val = suppReadMax + 1
                else :
                    max_val = suppRead[start_th_idx]
                min_val = suppRead[end_th_idx]

                elem = [x[1:]+[x[0]] for x in vc_result[sv_type] if (x[0] < max_val and x[0] >= min_val)]
                acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], elem, overlap, min_val) + [min_val] )
        else : 
            acc[sv_type] = {} 
    return acc

def getAccuracy(p, vc_result, overlap, minLen, maxLen, separation) :
    if separation == "b" :
        return getAccuracyByBin(p, vc_result, overlap, minLen, maxLen);        
    elif separation == "q" :
        return getAccuracyByQuantile(p, vc_result, overlap, minLen, maxLen);
    elif separation == "qc" :
        return getAccuracyByQuantileCum(p, vc_result, overlap, minLen, maxLen);
    else : # separation == "s" :
        return getAccuracyByScore(p, vc_result, overlap, minLen, maxLen);

def getAccuracyForCE(p, vc_result, overlap, minLen, maxLen, separation) :
    acc = {}
    for sv_type in [DUP_IMPR_STR, DEL_STR, INV_STR, BREAKPOINT_STR] : # for each SV type
        if ((sv_type in vc_result) == True) and len(vc_result[sv_type]) > 0 :

            acc[sv_type] = []

            # v_seqan contains [supp, se, pe, ce, re]
            for e_type in range(0,4) :

                # set evidence type
                _vc_result = copy.deepcopy(vc_result[sv_type])
                for j in range(0, len(_vc_result)) :
                    _vc_result[j][0] = vc_result[sv_type][j][0][e_type]

                suppRead = [ x[0] for x in _vc_result if x[0] > 0] # x[0] : num of supporting reads
                if len(suppRead) > 0 :
                    suppReadMin = min(suppRead)
                    suppReadMax = min( max(suppRead), MAX_SUPP )
                    elem = [x[1:]+[x[0]] for x in _vc_result if x[0] > 0]
                    acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], elem, overlap, 0)  + [e_type])
                else :
                    acc[sv_type].append( getAccuracyBySuppReads(sv_type, p[sv_type], [], overlap, 0)  + [e_type])
                    

            # uniq count
            se_only = 0
            pe_only = 0
            ce_only = 0
            se = 0
            pe = 0
            ce = 0
            for i in range(0, len(vc_result[sv_type])) :
                if vc_result[sv_type][i][0][1] > 0 and vc_result[sv_type][i][0][2] == 0 and vc_result[sv_type][i][0][3] == 0 :
                    se_only += 1
                if vc_result[sv_type][i][0][2] > 0 and vc_result[sv_type][i][0][1] == 0 and vc_result[sv_type][i][0][3] == 0 :
                    pe_only += 1
                if vc_result[sv_type][i][0][3] > 0 and vc_result[sv_type][i][0][1] == 0 and vc_result[sv_type][i][0][2] == 0 :
                    ce_only += 1
                se += vc_result[sv_type][i][0][1]
                pe += vc_result[sv_type][i][0][2]
                ce += vc_result[sv_type][i][0][3]

            # entire set
            acc[sv_type].append( [se, pe, ce, 0, 0, 0, 0, 0] )
            acc[sv_type].append( [se_only, pe_only, ce_only, 0, 0, 0, 0, 0] )
        

        else : 
            acc[sv_type] = {} 
    return acc


if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    sim_names = "|".join([key for key in c.CONFIG["genome-simulators"]])
    parser.add_argument("-b","--table", help="precision-recall(pr), soft-clipped evidences figure(ce), voting(vt)", default="pr", required=False)
    parser.add_argument("-t","--positive-type", help="[%s]" % sim_names, required=True)
    parser.add_argument("-f","--positive-file", help="A file that contains positive cases. This file should be produced by POSITIVE_TYPE", required=False)
    parser.add_argument("-p","--positive-file-parameters", help="Parameters for parsing positve files", default="", required=False)
    parser.add_argument("-g","--ref-genome-dir", help="A reference genome (with .fa). This must be exist in %s" % c.CONFIG['result-dir'], required=True)
    parser.add_argument("-o","--min-overlap", help="minimum reciprocal overlap [0 ~ 1.0]", required=True)
    parser.add_argument("-m","--min-sv-length", help="Minimum SV length in bp", default=-1, required=False)
    parser.add_argument("-a","--max-sv-length", help="Maxmimum SV length in bp", default=99999999999999999, required=False)
    parser.add_argument("-x","--max-threshold", help="Maxmimum threshold", default=100, required=False)
    parser.add_argument("-u","--min-supp", help="Minimum supp.", default=0, required=False)
    parser.add_argument("-j","--adjacency", help="maximum difference for adjacency", default=50, required=False)
    parser.add_argument("-s","--separation", help="separate by score(s), quantile(q), and bins(b)", default="s", required=False)
    args = parser.parse_args()

    # check dirs
    checkDirs()

    # params
    genomeDir = args.ref_genome_dir
    overlap = float(args.min_overlap)
    minLen = int(args.min_sv_length)
    maxLen = int(args.max_sv_length)
    maxTH = int(args.max_threshold)
    p_type = args.positive_type
    p_file = args.positive_file
    p_param = args.positive_file_parameters
    adjacency = args.adjacency
    separation = args.separation
    table = args.table
    minSupp = args.min_supp

    # open logger
    sys.stdout = Logger(c.CONFIG["log-dir"] + "/" + time.strftime("report_%Y%m%d.txt"))
    print("[START] " + " ".join(sys.argv))

    # load true positives
    if (p_type in dir()) == False :
        print("[ERROR] Can't find the function call : %s" % rm_name)
        exit(1)
    else :
        print("[FUNC] START: %s" % p_type)          
        if p_param == "" :
            p = run(p_type, [p_file] + [genomeDir])
        else :
            p = run(p_type, [p_file] + [p_param])
        print("[FUNC] END: %s" % p_type)

    # validation
    reportName = ""
    for param in sys.argv[1:]:
        if param[0] == "-" :
            reportName += param[1:] + ":"
        else :
            reportName += param + "#"
    reportName = reportName[:-1]
    print("[REPORT] : %s" % reportName)

    if table == "pr" :
        print("[VALIDATION] START: %s" % genomeDir)
        validation(p, genomeDir, overlap, adjacency, minLen, maxLen, maxTH, reportName, table, separation, minSupp)
        print("[VALIDATION] END: %s" % genomeDir)
    elif table == "ce" :
        print("[CLIPPED EVIDENCE] START: %s" % genomeDir)
        validation(p, genomeDir, overlap, adjacency, minLen, maxLen, maxTH, reportName, table, separation, minSupp)
        print("[CLIPPED EVIDENCE] END: %s" % genomeDir)
    elif table == "vt" :
        print("[VOTING] START: %s" % genomeDir)
        validation(p, genomeDir, overlap, adjacency, minLen, maxLen, maxTH, reportName, table, separation, minSupp)
        print("[VOTING END: %s" % genomeDir)

    print("[END]")