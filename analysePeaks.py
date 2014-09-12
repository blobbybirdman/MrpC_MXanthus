import os, sys, math, random
from operator import attrgetter

class AnalysePeaks:
    def __init__(self, input1, input2, outputF, ptt=None, cog=None):
        """Output root file name"""
        self.outputFile = output
        """Myxo genome size"""
        self.myxoGenomeSize = 9139763

        """Set thresholds"""
        self.replicateDistCutOff = 65
        self.upstreamCutOff = 400
        self.downstreamCutOff = 100

        """Create objects for gene and functional annotations"""
        self.ptt, self.cog = None, None
        if ptt:
            self.ptt = self.parseAnnotations(ptt)
        if cog:
            #self.cog = self.parseCOGAnnotations(cog)
            self.cog = self.parseHiggsAnnotations(cog)

        """Set up output file for summary information"""
        outputFile = open("".join((outputF, "-summary.txt")), 'w')
        
        """Read ChIP-seq replicate one and parse the file"""
        self.inputFile1 = input1
        outputFile.write("Reading input file: " + self.inputFile1 + '\n')
        self.peaks1 = self.parsePeaksFile(input1)
        outputFile.write("Peaks found: " + `len(self.peaks1)` + '\n')

        """Read Chip-seq replicate two and parse the file"""
        self.inputFile2 = input2
        outputFile.write("Reading input file: " + self.inputFile2 + '\n')
        self.peaks2 = self.parsePeaksFile(input2)
        outputFile.write("Peaks found: " + `len(self.peaks2)` + '\n')

        """Compare replicate datasets and identify peaks colocated with no cutoff distance"""
        outputFile.write("Searching for replicate peaks with no distance cutoff\n")
        self.unfilteredReplicatePeaks = self.findSharedPeaks(self.peaks1, self.peaks2)
        outputFile.write("Number of shared peaks identified: " + `len(self.unfilteredReplicatePeaks)` + '\n\n')

        """Write distances between unfiltered replicate peaks to output file"""
        self.writePeakDistances(self.outputFile + "-unfilteredPeakDistances-noAnnot.txt", self.unfilteredReplicatePeaks)

        """Compare replicate datasets and identify peaks colocatedwithin cutoff distance"""
        outputFile.write("Searching for peaks colocated within " + `self.replicateDistCutOff` + "bp across replicates\n")
        self.replicatePeaks = self.findSharedPeaks(self.peaks1, self.peaks2, self.replicateDistCutOff)
        outputFile.write("Number of shared peaks identified: " + `len(self.replicatePeaks)` + '\n\n')

        """Write distances between replicate peaks to output file"""
        self.writePeakDistances(self.outputFile + "-peakDistances-noAnnot.txt", self.replicatePeaks)

        """Generate a set of peaks with randomised locations in myxo genome. Excess needed to ensure sufficient after filtering"""
        randomPeaks = self.getRandomPeaks(1608)
        
        """Copy replicatePeaks and generate matching randomPeaks so can explore annotations without filtering (without refactoring completely *workaround*)"""
        self.origReplicatePeaks = self.replicatePeaks[:]
        origRandomPeaks = self.getRandomPeaks(len(self.origReplicatePeaks))

        """Convert random peaks to replicate peak format for compatibility"""
        self.randomPeaks = self.findSharedPeaks(randomPeaks, randomPeaks, self.replicateDistCutOff)
        self.origRandomPeaks = self.findSharedPeaks(origRandomPeaks, origRandomPeaks, self.replicateDistCutOff)

        """Annotate peaks (random and data) without filtering"""
        self.origCoding = len(filter(self.isCoding, self.origReplicatePeaks))
        print "Replicate peaks:"
        print "Coding: ", self.origCoding
        print "Non-coding: ", len(self.origReplicatePeaks) - self.origCoding

        self.ranCoding = len(filter(self.isCoding, self.origRandomPeaks))
        print "Random"
        print "Coding: ", self.ranCoding
        print "Non-coding: ", len(self.origRandomPeaks) - self.ranCoding

        outputFile.write(`self.origCoding` + " unfiltered replicate peaks found in coding regions. " + `len(self.origReplicatePeaks)-self.origCoding` + " found in non-coding regions.\n")
        outputFile.write(`self.ranCoding` + " unfiltered random peaks found in coding regions. " + `len(self.origRandomPeaks)-self.ranCoding` + " found in non-coding regions.\n\n")


        """Annotate each peak with annotations of nearby ORFs (default -400 to +100)"""
        self.replicatePeaks = self.annotatePeaks(self.replicatePeaks, self.upstreamCutOff, self.downstreamCutOff)
        outputFile.write(`len(self.replicatePeaks)` + " filtered replicate peaks found within " + `self.upstreamCutOff`)
        outputFile.write("bp upstream and " + `self.downstreamCutOff` + "bp downstream of ORF annotations\n")

        """Calculate number of replicate peaks located between divergent promoters"""
        numAnnotations = [len(peak.annotations) for peak in self.replicatePeaks if len(peak.annotations) > 1 and len(peak.annotations) < 3]
        outputFile.write(`len(numAnnotations)` + " of these filtered peaks are located in divergent promoters with ambiguous association\n\n")
        
        """Write details of paired peaks from replicates"""
        self.writePeakPairs(self.outputFile, self.replicatePeaks)

        """Write frequency distribution data for replicate paired peaks"""
        self.writePairDistanceDistribution(self.outputFile, self.replicatePeaks)

        
        
        """Annotate randomized peaks and then select a number equal to the number of replicate peaks after filtering steps"""
        self.randomPeaks = self.annotatePeaks(self.randomPeaks, self.upstreamCutOff, self.downstreamCutOff)
        
        outputFile.write(`len(self.randomPeaks)` + " filtered randomized peaks found within " + `self.upstreamCutOff`)
        outputFile.write("bp upstream and " + `self.downstreamCutOff` + "bp downstream of ORF annotations\n")

        """Calculate number of random peaks located between divergent promoters"""
        numAnnotations = [len(peak.annotations) for peak in self.randomPeaks if len(peak.annotations) > 1]
        outputFile.write(`len(numAnnotations)` + " of these filtered peaks are located in divergent promoters with ambiguous association\n\n")

        """Write out distances from replicate peaks to nearest ORFs to output file"""
        self.writeORFDistances(self.outputFile + "-peakORFDistances.txt", self.replicatePeaks)

        """Write out distance distribution from replicate peaks to nearest ORFs to output file"""
        self.writeORFDistanceDistribution(self.outputFile + "-peakORFDistanceDistribution.csv", self.replicatePeaks, \
                                              self.upstreamCutOff, self.downstreamCutOff)

        """Write out distances from randomized peaks to nearest ORFs to output file"""
        self.writeORFDistances(self.outputFile + "-randomORFDistances.txt", self.randomPeaks)

        """Write out distance distribution from randomized peaks to nearest ORFs to output file"""
        self.writeORFDistanceDistribution(self.outputFile + "-randomORFDistanceDistribution.csv", self.randomPeaks, \
                                              self.upstreamCutOff, self.downstreamCutOff)
        
        """Run r script to generate plots of various distributions"""
        os.system("./generatePeakORFStats.r &>/dev/null")


        self.intergenicDistribution = self.getPeakGenomicDistribution(self.replicatePeaks)
        outputFile.write("Genomic analysis of orig peak distribution:\n")
        outputFile.write('upstream:' + `len(self.intergenicDistribution['IU'])` + '\n')
        outputFile.write('coding:' + `len(self.intergenicDistribution['CO'])` + '\n')

        self.randomIntergenicDistribution = self.getPeakGenomicDistribution(self.randomPeaks)
        outputFile.write("Genomic analysis of orig randomized peak distribution:\n")
        outputFile.write('upstream: ' + `len(self.randomIntergenicDistribution['IU'])` + '\n')
        outputFile.write('coding: ' + `len(self.randomIntergenicDistribution['CO'])` + '\n\n')

        outputFile.write("\n\npeaks: " + `len(self.replicatePeaks)`)
        outputFile.write("\tannotations: " + `len(self.intergenicDistribution['IU']) + len(self.intergenicDistribution['CO'])`)
        outputFile.write("\nRandom peaks: " + `len(self.randomPeaks)`)
        outputFile.write("\tannotations: " + `len(self.randomIntergenicDistribution['IU']) + len(self.randomIntergenicDistribution['CO'])`)

        """Write summary of distribution of peaks between ORFs and intergenic regions"""
        self.writePeakDistribution(self.outputFile + "-intergenicDistribution", self.intergenicDistribution)

        """Write summary of distribution of peaks between ORFs and intergenic regions"""
        self.writePeakDistribution(self.outputFile + "-randomIntergenicDistribution", self.randomIntergenicDistribution)

        """Write COG annotation summary"""
        self.writeCOG(output)

        """Write MXAN numbers"""
        self.writeMXAN(output)

        """Write output"""
        self.writePeaks(output)


    def parsePeaksFile(self, f):
        """Iteration through QuEST output file extracting individual peaks"""
        peaks = []
        file = open(f, 'r')
        for line in file.readlines():
            if len(line) == 1 or line[0] == 'R':
                """Ignore lines that do not contain information pertaining to individual peaks"""
                pass
            else:
                peaks.append(Peak(line.split()))

        """sort peaks by coordinate"""
        peaks = sorted(peaks, key=attrgetter('coord'))
        return peaks


    def parseAnnotations(self, ptt):
        """Iterate through genome annotation file generating an object for each gene annotation"""
        file = open(ptt, 'r')
        annotations = []
        for line in file.readlines()[3:]:
            annotations.append(PTT(line))
        return annotations

   
    def parseCOGAnnotations(self, cog):
        """Iterate through functional annotation file generating an object for each function annotation"""
        file = open(cog, 'r')
        cogAnnotations = {}
        for line in file.readlines():
            items = line.split()
            if len(items[-1]) > 1:
                """ensure COG annotation is available for ORF"""
                MXAN = items[-2].split('_')[1]
                cog = items[-1][-1]
                cogAnnotations[MXAN] = COGAnnotation(MXAN, items[-2], cog)
        return cogAnnotations

    def parseHiggsAnnotations(self, cog):
        """Iterate through Higg's functional annotation file generating an object for each function annotation"""
        file = open(cog, 'r')
        cogAnnotations = {}
        for line in file.readlines()[1:]:
            items = line.split('\t')
            MXAN = items[0][4:]
            gene = items[2]
            product = items[1]
            cog = items[3]
            cogSubCat = items[4]
            regulation = items[5]
            cogAnnotations[MXAN] = COGAnnotation(MXAN, gene, cog)
        return cogAnnotations

    
    def findSharedPeaks(self, p1, p2, cutoff=None):
        """Find all nearest peaks, and all nearest peaks within a predetermined cutoff distance"""
        replicatePeaks = []
        for peak1 in p1:
            peaks = []
            bestDistance = None
            bestPeak = None
            for peak2 in p2:
                """Calculate distance between current peaks"""
                distance = peak1.coord-peak2.coord

                """Check if current pairing is the closest found yet"""
                if bestDistance == None or math.fabs(distance) < math.fabs(bestDistance):
                    bestDistance = distance
                    bestPeak = peak2

            """If closest peak pair found is under cutoff then include in closePeaks"""
            if cutoff:
                if math.fabs(bestDistance) < cutoff:
                    replicatePeaks.append(ReplicatePeaks([peak1, bestPeak], bestDistance))
            else:
                """No cutoff used primarily during testing"""
                replicatePeaks.append(ReplicatePeaks([peak1, bestPeak], bestDistance))
        
        """sort replicate peaks by average peak score"""
        replicatePeaks = sorted(replicatePeaks, key=attrgetter('avgScore'))
        replicatePeaks.reverse()

        return replicatePeaks


    def annotatePeaks(self, peaks, upstreamCutoff=None, downstreamCutoff=None):
        filteredPeaks = []
        for peak in peaks:
            coord = peak.coord 
            distance = None

            """Iterate through gene annotations and find those within range for current peak"""
            for ptt in self.ptt:
                """Calculate distance from peak to start of current annotation"""
                distance = math.fabs(ptt.start - coord)

                """check if cog annotation available"""
                cog = None
                try:
                    cog = self.cog[ptt.synonym]
                except:
                    pass
                
                """If annotation is upstream"""
                if (ptt.strand == '+' and ptt.start >= coord) or \
                        (ptt.strand == '-' and ptt.start <= coord):
                    """If annotation is within upstream cutoff add to annotations"""
                    if distance <= upstreamCutoff:
                        peak.annotations.append(Annotation(ptt, '-', distance, cog))
                else:
                    """if annotation is within downstream cutoff"""
                    if distance <= downstreamCutoff:
                        peak.annotations.append(Annotation(ptt, '+', distance, cog))
            if len(peak.annotations) > 2:
                newAnnot = {}
                for annotation in peak.annotations:
                    newAnnot[annotation.distance] = annotation
                peak.annotations = []
                keys = newAnnot.keys()
                keys.sort()
                for each in keys[:2]:
                    peak.annotations.append(newAnnot[each])
            if peak.annotations:
                filteredPeaks.append(peak)
            
        return filteredPeaks

    def isCoding(self, peak):
        coord = peak.coord 

        """Iterate through gene annotations and check if peak falls in a coding region"""
        for ptt in self.ptt:
            if ptt.strand == '+':
                if coord >= ptt.start and coord <= ptt.end:
                    return True
            else:
                if coord >= ptt.end and coord <= ptt.start:
                    return True
        return False

    def writePeaks(self, output):
        """Generate output file showing annotation details for each significant ChIP-seq peak"""
        file = open("".join((output, "-peaks.csv")), 'w')
        file.write("Rank\tCoord\tPair Distance\t# annotations matching\tMXAN")
        file.write("\tGene Product\tORF Start\tDistance to ORF Start\tPeak Postion Relative to ORF Direction\n")
        rank = 0
        for pair in self.replicatePeaks:
            rank += 1
            if len(pair.annotations) == 0:
                """if no annotations associated with peak pair then just print pair summary"""
                file.write(`rank` + '\t' + `pair.coord` + '\t' + `abs(pair.distance)`)
                file.write('\t' + `len(pair.annotations)` + '\n')
            for annotation in pair.annotations:
                if len(pair.annotations) > 2:
                    continue
                """if annotations associated with pair, print summary line for each annotation"""
                file.write(`rank` + '\t' + `pair.coord` + '\t' + `abs(pair.distance)`)
                file.write('\t' + `len(pair.annotations)`)
                file.write('\t' + annotation.ptt.synonym + '\t' + annotation.ptt.product + '\t' + `annotation.ptt.start
                           ` + '\t' + `int(annotation.distance)` + '\t' + annotation.direction)
                if annotation.cog:
                    """If COG category for annotation also know, include that in summary line"""
                    file.write('\t' + annotation.cog.cogCat)
                file.write('\n')
        file.close()

    def writeMXAN(self, output):
        file = open("".join((output, "-MXAN.txt")), 'w')
        for peak in self.replicatePeaks:
            for annotation in peak.annotations:
                file.write(annotation.ptt.synonym)
                file.write('\n')
        file.close()

    def writePeakPairs(self, output, peakPairs):
        """Generate output file showing details for both records in matched replicate pairs"""
        file = open("".join((output, "-pairs.txt")), 'w')
        for pairs in peakPairs:
            file.write('start\n')
            file.write(`pairs.peaks[0]`+'\n')
            file.write(`pairs.peaks[1]`+'\n')
            file.write("Distance: " + `int(math.fabs(pairs.distance))` +'\n')
        file.close()

    def writePairDistanceDistribution(self, output, pairs):
        """Generate output file with the distance distribution between replicates for all peaks and thresholded peaks"""
        file = open("".join((output, "-distance.csv")), 'w')
        freq = {}
        distances = []
        for peak in pairs:
            """Append distance between replicate peaks to list"""
            distance = peak.distance
            distances.append(distance)
            """Increment dictionary of frequencies with new distance"""
            freq[distance] = 1 + freq.get(distance, 0)
        for item in sorted(freq.items()):
            """Write line for each distance with a frequency of one or more"""
            file.write(`item[0]` + '\t' + `item[1]` + '\n')
        file.close()

        """Generate output file summarizing the average distance  and stdDev between replicate peaks"""
        file = open("".join((output, "-distance.txt")), 'w')
        mean, stdDev = self.getStdDev(distances)
        file.write("Mean distance: " + `mean` + '\n')
        file.write("StdDev: " + `stdDev` + '\n')
        file.close()

    def getStdDev(self, items):
        """calculate average distance between peaks"""
        average = sum(items)*1.0/len(items)
        """Calculate the variance of the average distance between peaks"""
        varience = map(lambda x: (x-average)**2, items)
        """Calculate the stdDev from the variance"""
        avgVar = sum(varience)/len(varience)
        stdDev = math.sqrt(avgVar)
        return average, stdDev

    def getRandomPeaks(self, numPeaks):
        """Generate a number of random peaks, picking a random coordinate in the myxo geneome for each"""
        peaks = []
        for x in range(numPeaks):
            """Generate random genome coordinate"""
            coord = random.randint(1, self.myxoGenomeSize)
            """Generate Peak object using random genome coordinate as a basis"""
            peaks.append(Peak(["random", "M_xanthus", coord, 0, 0, 0, 0, 0, `coord`+'-'+`coord`,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        return peaks

    def writeORFDistances(self, output, peaks):
        """Generate output file listing distances between peaks and any associated annotation, one per line"""
        file = open(output, 'w')
        for peak in peaks:
            if peak.annotations:
                for annotation in [peak.annotations[0]]:
                    """Report distance for all annotations associated with each peak"""
                    distance = int(annotation.distance)
                    if (annotation.direction == '-'):
                        """If annotation is downstream of ORF start, ensure distance is reported as negative"""
                        distance = -distance
                    file.write(`distance` + '\n')
        file.close()

    def writePeakDistances(self, output, peaks):
        """Generate output file listing distances between peaks, one per line"""
        file = open(output, 'w')
        for peak in peaks:
                file.write(`int(peak.distance)` + '\n')
        file.close()

    def writeORFDistanceDistribution(self, outputFile, peaks, upstreamCutOff, downstreamCutOff):
        """Calculate frequency bins and write distribution of peak to ORF start distribution to file"""
        distribution = {}
        file = open(outputFile, 'w')
        interval = ((upstreamCutOff+downstreamCutOff)/50)
        bins = [x for x in range(-1*upstreamCutOff, downstreamCutOff, interval)]
        for peak in peaks:
            if peak.annotations:
                for annotation in [peak.annotations[0]]:
                    distance = int(annotation.distance)
                    direction = annotation.direction
                    if (direction == '-'):
                        distance = -distance
                    for eachBin in bins:
                        if distance < (eachBin + interval -1):
                            distribution[eachBin] = 1 + distribution.get(eachBin, 0)
                            break
        keys = distribution.keys()
        keys.sort()
        for key in keys:
            file.write(`key` +  "\t" + `distribution[key]` + '\n')            
        file.close()

    def writeCOG(self, output):
        """Generate output file summarising cog categories for genes associated with ChIP-seq peaks"""
        file = open("".join((output, "-cog.txt")), 'w')
        genomeCogCatTable = {}
        cogCatTable = {}
        cogSum = 0.
        for pair in self.replicatePeaks:
            for annotation in pair.annotations:
                """For each annotation associated with each peak count frequency of COG categories represented"""
                if annotation.cog:
                    cogSum += 1
                    cogCatTable[annotation.cog.cogCat] = 1 + cogCatTable.get(annotation.cog.cogCat, 0)
        
        genomeCogSum = 0.
        for key in self.cog.keys():
            genomeCogSum += 1
            genomeCogCatTable[self.cog[key].cogCat] = 1 + genomeCogCatTable.get(self.cog[key].cogCat, 0)
            
        file.write("CogCat\tGenome\tChIP\tOther\tExpected\tExpected other\n")
        for key in genomeCogCatTable.keys():
            try:
                cogCat = cogCatTable[key]
                other = int(cogSum-cogCat)
            except:
                cogCat = 0
                other = cogSum
            genomeCat = genomeCogCatTable[key]
            expected = genomeCat*cogSum/genomeCogSum
            expectOther = cogSum-expected

            file.write(key + '\t' + `genomeCat` + '\t')
            file.write(`cogCat` + '\t' + `other` + '\t' + `expected` + '\t' + `expectOther` + '\n')
        file.write('total:\t' + `sum(genomeCogCatTable.values())` + '\t' + `sum(cogCatTable.values())` + ' \n')
        file.close()

        

    def getPeakGenomicDistribution(self, peaks):
        """Determine how many peaks fall in either coding regions, or upstream/downstream intergenic regions"""
        distribution = {"CO": [], "IU": [], "ID": []}
        for peak in peaks:
            if len(peak.annotations) > 2:
                continue
            for annotation in peak.annotations:
                start, end = annotation.ptt.start, annotation.ptt.end
                if annotation.ptt.strand == '-':
                    start, end = end, start
                if annotation.direction == '-':
                    distribution['IU'].append(int(annotation.distance))
                elif annotation.direction == '+' and (start <= peak.coord <= end):
                    distribution['CO'].append(int(annotation.distance))
                else:
                    distribution['ID'].append(int(annotation.distance))
        return distribution

    def writePeakDistribution(self, output, distribution):
        """Write summary of distribution of peaks among coding and non-coding regions"""
        file = open("".join((output, '-summary.txt')), 'w')
        file.write("Summary\n=======\n")
        file.write("Peaks in Upstream intergenic regions: " + `len(distribution['IU'])` + "\n")
        file.write("Peaks in ORFs: " + `len(distribution['CO'])` + "\n")
        file.write("Peaks in Downstreams intergenic regions: " + `len(distribution['ID'])` + "\n\n\n")
        file.close()

 
class ReplicatePeaks:
    """Class for storing peaks found in both experimental replicates"""
    def __init__(self, peaks, distance):
        """Initialise variables"""
        self.peaks = peaks
        self.distance = int(distance)
        self.coord = self.getCoord()
        self.annotations = []
        self.avgScore = self.getAvgScore()
        self.avgEF = self.getAvgEF()

    def getCoord(self):
        """Find average coord for the replicated peaks"""
        total = 0
        for peak in self.peaks:
            total += peak.coord
        return int(total/len(self.peaks))

    def getAvgScore(self):
        """Find average score associated with replicate peaks"""
        total = 0
        for peak in self.peaks:
            total += peak.qv
        return float(total/len(self.peaks))

    def getAvgEF(self):
        """Find average enrichment (relative to control) of replicate peaks"""
        sum = 0
        for peak in self.peaks:
            sum += peak.ef
        return float(sum/len(self.peaks))

    def __repr__(self):
        """representation of class for printing"""
        return "<%s %s %s %s>" % (self.peaks[0].id, self.peaks[1].id, self.coord, self.annotations)


class Annotation:
    """Class for holding peak annotation data"""
    def __init__(self, annotation, direction, distance, cog):
        self.ptt = annotation
        self.direction = direction
        self.distance = distance
        self.cog = cog

    def __repr__(self):
        """representation of class for printing"""
        return "<%s %s %i, %s>" % (self.ptt, self.direction, self.distance, self.cog)


class PTT:
    """Class for holding Gene annotation data"""
    def __init__(self, record):
        """Extract data from record"""
        self.location, self.strand, self.length, self.PID, self.gene, self.synonym, self.code, self.COG, self.product =  record.split('\t')
        
        """"Trim newline from product string"""
        self.product = self.product[:-1]

        """extract MXAN number from synonym string"""
        self.synonym = self.synonym.split('_')[1]

        """Extract (int) start and end values from location, ensuring that start is smallest coordinate on + strand"""
        start, end = self.location.split('..')
        if self.strand == '+':
            self.start, self.end = int(start), int(end)
        elif self.strand == '-':
            self.start, self.end = int(end), int(start)
        else:
            print "improper strand"

    def __repr__(self):
        """representation of class for printing"""
        return "<%s %i %s %s>" % (self.synonym, self.start, self.strand, self.gene)


class COGAnnotation:
    """Class for holding functional annotation data"""
    def __init__(self, gene, cogRef, cogCat):
        self.gene, self.cogRef, self.cogCat = gene, cogRef, cogCat

    def __repr__(self):
        """representation of class for printing"""
        return "<%s %s %s>" % (self.gene, self.cogRef, self.cogCat)

class Peak:
    """Basic class for holding information pertaining to an individual peak"""
    def __init__(self, record):
        """Read in split() line from Quest output and allocate to class variables"""
        self.id, self.species, self.coord, self.count, self.controlCount, self.range, \
            self.ef, self.ps, self.cor, self.qv, self.pv, self.qvRank = record[0], record[1], \
            int(record[2]), float(record[4]), float(record[6]), record[8], float(record[10]), float(record[12]), \
            float(record[14]), float(record[16]), float(record[18]), float(record[20])

        """extract stard and end coords from range field"""
        self.rangeStart, self.rangeEnd = self.range.split('-')
    
    def __repr__(self):
        """representation of class for printing"""
        return "<%s %s %s %s %s %s %s %s %s %s %s %s>" % (self.id, self.species, self.coord, self.count, \
                                                          self.controlCount, self.range, self.ef, self.ps, \
                                                              self.cor, self.qv, self.pv, self.qvRank)


if __name__ == '__main__':
    ptt = None
    cog = None
    if len(sys.argv) == 5:
        """Input, output, and gene annotation files"""
        f1, f2, output, ptt = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    elif len(sys.argv) == 4: 
        """Input and output files"""
        f1, f2, output = sys.argv[1], sys.argv[2], sys.argv[3]
    elif len(sys.argv) == 6:
        """Input, output, gene annotationa and functional category files"""
        f1, f2, output, ptt, cog = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    else:
        print "Num of Args: ", len(sys.argv)
        print "Error: A file containing paired replicate QuEST peak calls are required to be supplied on the command line for comparison along with an output file"
        sys.exit()

    AnalysePeaks(f1, f2, output, ptt, cog)
