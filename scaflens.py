#!/usr/bin/env python

# output scaffold names and lengths
# 17Dec14 JBH adding several new options:
#    -t or -n to output, instead of each scaffold name/len, either total of lens or with -n total and the scaff N50 cal
#    either -t or -n can be followed by a list of numbers (eg 1000 5000 10000), this will cause additional
#    total lens and N50 values to be output where scaffolds of each length in args are excluded from the calculation
# 01Jan15 JBH add -r option to show ranges of sizes of the scaffolds and ranges of sizes of the singletons

import sys, os, re
from collections import Counter

if len(sys.argv) < 2:
    print "\nUsage: scaflens.py <scafSeq_file_from_SOAP_assembly> [ [-t|-n [-C] [scaff_len ...]] | -r ]\n"
    print "         Outputs scaffold names and their sizes in a list (unless -t, -n or -r option used)."
    print "         You can pipe output to sort -nr -k2,2 to see the scaffolds listed largest to smallest.\n"
    print "         -t option shows total length of scaffolds (includes Ns). Unscaffolded singleton contigs are excluded (unless -C used)."
    print "         -n shows the N50 count along with the total."
    print "         To see N50 and/or total count excluding scaffolds of a particular size or shorter"
    print "            you can add a number (or list of numbers) after the -t or -n option"
    print "         -C without this only scaffolds are included, using this also includes singleton contigs;"
    print "            when using -C 200 contigs/scaffolds 200 or smaller are ignored\n"
    print "         -r shows ranges of scaffold sizes and ranges of singleton contig sizes (other options ignored)\n"
    print "\n         Examples: scaflens.py soap_asm12.scafSefSeq"
    print "                   scaflens.py soap_asm12.fasta | sort -nr -k2,2"
    print "                   scaflens.py soap_asm12.fasta -n50 9999 99999"
    print "                   scaflens.py *.fasta -n -C 500 1000 10000 100000 1000000"
    print "                   scaflens.py -r *.fasta\n"
    sys.exit(0)
    
def showInfo(excludeSize, usingContigs=False):
    global scafSizeList # this is assumed to be sorted largest to smallest sizes 
    global len_singletons_le_200, num_singletons_le_200 # if useContigs these will lump a lot together
    
    totalScafs = len(scafSizeList)
    lngth = 0; numScafs = 0
    for s in xrange(totalScafs):
        sz = scafSizeList[s]
        if excludeSize > 0 and sz <= excludeSize: # excludeSize of 0 means skip test, so we can count all
            break
        numScafs +=  1
        lngth += scafSizeList[s]
        
    # only in case of usingContigs==True and counting all (ie excludeSize == 0), we add back the size and number of
    # singleton contigs 200 of less in size
    if usingContigs and excludeSize==0:
        numScafs += num_singletons_le_200
        lngth += len_singletons_le_200
    
    scafstr = "scaffolds"
    excludeStr = "(excludes singletons)"
    if usingContigs:
        scafstr = "scaffolds/contigs"
        excludeStr = "(includes singleton Contigs)"
    if excludeSize > 0:
        excludeStr = "(excluding "+scafstr+" " + str(excludeSize) + " or shorter)"
        
    if not showN50: 
        print "Length", lngth, excludeStr
        return
    
    N50 = 0
    midPoint = lngth / 2
    BPsofar = 0
    for s in xrange(totalScafs): # looping from biggest scaffold lengths to smallest, see when we cross 50%
        BPsofar += scafSizeList[s]
        if BPsofar > midPoint:
            N50 = scafSizeList[s] # this is the length of the contig that sent us over the halfway point
            break
    
    print "N50", N50, "Length", lngth, "in " + str(numScafs), scafstr, excludeStr
#end showInfo()

quartiles100 = [0,0,0,0] # track the sizes 100-124, 125-149, 150-174, 175-199
quartiles200 = [0,0,0,0]
def getRanges(scafFile):
    def sizeToRangeInt(sz):
        if sz > 499: # just use the size for the counter (not that many 500 or over)
            return sz
        if sz >= 400: # 400-499 range use 400
            return 400
        if sz >= 300: # 300-399 range use 300
            return 300
        if sz >= 200: # 200-299 range use 200
            qix = (sz-200) / 25
            quartiles200[qix] += 1
            return 200
        if sz >= 100: # 100-199 range use 100
            qix = (sz-100) / 25
            quartiles100[qix] += 1
            return 100
        if sz >= 50: # 50-99 range use 50
            return 50
        else: # 1-49 range use 25
            return 25
        
    fh = open(scafFile)
    
    scafSizes = Counter()
    singletonSizes = Counter()
    singletonContig = True # set this to false when we see an N in a line
    size = 0
    ln = fh.readline().strip()
    while ln:
        if ln[0]=='>':
            if size > 0:
                if singletonContig:
                    rngint = sizeToRangeInt(size)
                    singletonSizes[rngint] += 1
                else:
                    scafSizes[size] += 1
            size = 0
            singletonContig = True
        else:
            if singletonContig and ln.find("N") != -1:
                singletonContig = False
            size += len(ln)
        ln = fh.readline().strip()
    fh.close()
    
    if size > 0:
        if singletonContig:
            rngint = sizeToRangeInt(size)
            singletonSizes[rngint] += 1
        else:
            scafSizes[size] += 1
            
    return scafSizes, singletonSizes
#end getranges()

def printRanges(scafSizes, singletonSizes):
    def intToRangeStr(rngint):
        if rngint < 500:
            if rngint >= 100:
                strint = str(rngint)
                return strint[0] + "00-" + strint[0] + "99  "
            elif strint >= 50:
                return "50-99"
            else:
                return "1-49"
        elif rngint < 1000:
            return "500-999  "
        elif rngint < 2000:
            return "1000-1999"
        elif rngint < 5000:
            return "2000-4999"
        elif rngint < 10000:
            return "5000-9999"
        elif rngint < 100000: #10,000 - 99,999, see which 10 tho it is
            tho = 20
            while tho <= 100:
                if rngint < (tho*1000):
                    return str(tho-10)+"000-" + str(tho-10)+"999"
                tho += 10
            return "10000-99999" # fail safe
        elif rngint < 1000000:
            return "100000-999999"
        else:
            return "1000000+"
    
    def printSizeCounter(sizeCounter, msg, showQuartile=False):
        range = Counter()
        for k,v in sizeCounter.iteritems():
            r = intToRangeStr(k)
            range[r] += v
        
        items = [(k, v) for k, v in range.items()]
        sortedRanges = sorted(items, key=lambda item: int(re.match(r'\d+', item[0]).group())) # return first portion of range string as an int for sorting
        maxCount = max(items, key=lambda x: x[1])[1]; maxCountDigits = len(str(maxCount))
        print "\n"+msg
        print "Length Ranges\tCounts"
        for itm in sortedRanges:
            count = str(itm[1])
            print itm[0] + '\t' + count.rjust(maxCountDigits),
            if showQuartile:
                if itm[0].strip() == "100-199":
                    print " " + quartileStr(quartiles100, 100),
                elif itm[0].strip() == "200-299":
                    print " " + quartileStr(quartiles200, 200),
            print
    
    def quartileStr(qlist, base):
        qs = "("
        for qix in xrange(4):
            qtrStart = base + qix*25
            qs += str(qtrStart) + "-" + str(qtrStart+24) + " " + str(qlist[qix])
            if qix < 3:
                qs += ", "
        return qs + ")"
    
    printSizeCounter(scafSizes, "Scaffold")
    printSizeCounter(singletonSizes, "Singleton Contig", True)
    
#end printranges()
    
scafFile = ""
al = re.compile('>[a-zA-Z]+')
scafSizeList = []
sortLens = False; showTotals = False; showN50 = False; showRanges = False; scaffLenCutOffs = []
includeContigs = False
len_singletons_le_200 = 0 # this var and next only used if includeContigs is True
num_singletons_le_200 = 0

numArgs = len(sys.argv)
ixArg = 0
while ixArg < (numArgs-1):
    ixArg += 1
    arg = sys.argv[ixArg]
    if arg[0] != "-":
        if scafFile == "":
            scafFile = arg
    else: # arg starts with hyphen
        if arg[1] == "t" or arg[1] == "n": # it's a -t or -n argument
            showTotals = True
            if arg[1] == "n":
                showN50 = True
            if scafFile != "": # already have file, look for list of numbers (or -C) after -n or after -t
                while ixArg < (numArgs-1):
                    ixArg += 1
                    arg = sys.argv[ixArg]
                    if arg.isdigit():
                        scaffLenCutOffs.append(int(arg))
                    elif arg[:2] == "-C":
                        includeContigs = True
        elif arg == "-r": # show Ranges
            showRanges = True

if showRanges:
    scafSizes, singletonSizes = getRanges(scafFile)
    printRanges(scafSizes, singletonSizes)
else:
    fh = open(scafFile)
    ln = fh.readline().strip()
    lastWasScaffold = False
    size = 0
    while ln:
        if ln[0]=='>':
            if lastWasScaffold:
                scafSizeList.append(size)
            elif includeContigs: # special handling if we are looking at contigs
                if size > 200: # just treat it as though it was a scaffold
                    scafSizeList.append(size)
                elif size > 0:
                    len_singletons_le_200 += size
                    num_singletons_le_200 += 1
            size = 0
            lastWasScaffold = False
            m = al.match(ln)
            if m:
                txt = m.group()[1:]
                lastWasScaffold = (txt == 'scaffold' or txt == 'Scaffold') # 24Dec14 JBH looks like SOAP1 used uppercase S
        else:
            size += len(ln)
        ln = fh.readline().strip()
    fh.close()

    if lastWasScaffold:    
        scafSizeList.append(size)
    elif includeContigs: # special handling if we are looking at contigs
        if size > 200: # just treat it as though it was a scaffold
            scafSizeList.append(size)
        elif size > 0:
            len_singletons_le_200 += size
            num_singletons_le_200 += 1
        
    if not showTotals: # common case, output list of scaffolds with their sizes
        for s in xrange(len(scafSizeList)):
            print "scaffold"+str(s+1), scafSizeList[s]
    else: # we're going to show totals and possibly N50s, for all scaffolds and potentially excluding some of the smaller ones
        scafSizeList.sort(reverse=True)
        totalLen = sum(scafSizeList)
        showInfo(0, includeContigs) # show info for all the scaffolds (but exclude the singletons)
        if len(scaffLenCutOffs) > 0:
            for ex in scaffLenCutOffs:
                if not includeContigs or ex > 200: # don't allow length exclusions 200 or less if using singleton contigs
                    showInfo(ex, includeContigs)
