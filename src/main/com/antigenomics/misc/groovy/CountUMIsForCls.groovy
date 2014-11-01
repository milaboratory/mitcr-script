package mitcr

@Grab(group = 'com.milaboratory', module = 'mitcr-groovy', version = '1.0.3')
import com.milaboratory.core.clone.Clone
import com.milaboratory.core.clone.CloneSetClustered
import com.milaboratory.core.segment.AlleleSet
import com.milaboratory.core.segment.Gene
import com.milaboratory.core.segment.SegmentSet
import com.milaboratory.core.segment.Species
import com.milaboratory.core.sequencing.io.fastq.SFastqReader
import com.milaboratory.mitcr.clonegenerator.SequencingReadLink
import com.milaboratory.mitcr.pipeline.DefaultAnalysisListener
import com.milaboratory.mitcr.pipeline.FullPipeline
import com.milaboratory.mitcr.pipeline.ParameterPresets
import com.milaboratory.mitcr.qualitystrategy.IlluminaQualityInterpretationStrategy
import com.milaboratory.util.CompressionType
import com.milaboratory.util.ProgressReporter


if (args.length < 2) {
    println "groovy CountUMIsForCls input_R2 output NNN_COUNT_THRESHOLD SPECIES GENE QUALITY"
    println "Allowed values: NNN_COUNT_THRESHOLD=1..inf SPECIES=HomoSapines,MusMusculus GENE=TRA,TRB QUALITY=2..40"
    println "Default values: NNN_COUNT_THRESHOLD=1 SPECIES=HomoSapines GENE=TRB QUALITY=20"
    System.exit(0)
}

def r2FileName = args[0], outFileName = args[1]
def nnnCountThreshold = Integer.parseInt(args.length > 2 ? args[2] : "1")
def species = args.length > 3 ? args[3] : "HomoSapiens", gene = args.length > 4 ? args[4] : "TRB"
def quality = Integer.parseInt(args.length > 5 ? args[5] : '20')

def nnnList = new ArrayList<String>()
def reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
def read
while ((read = reads.take())) {
    nnnList.add(read.getDescription().split(" UMI:")[1].split(":")[0])
}

def params = ParameterPresets.jPrimer
params.species = Species."$species"
params.gene = Gene."$gene"
params.qualityInterpretationStrategy = new IlluminaQualityInterpretationStrategy((byte) quality)

reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
def pipeline = new FullPipeline(reads, params, true)
def writeAlleles = { AlleleSet alleles ->
    alleles.collect { it.fullName }.join(",")
}
def writeSegments = { SegmentSet segments ->
    segments.collect { it.segmentName }.join(",")
}
pipeline.setAnalysisListener(new DefaultAnalysisListener() {
    @Override
    void afterClusterization(CloneSetClustered clusterizedCloneSet) {
        new File(outFileName).withPrintWriter { pw ->
            pw.println(["NNNs", "Count", "Part",
                    "CDR3 nucleotide sequence", "CDR3 nt seq quality",
                    "CDR3 amino acid sequence",
                    "V alleles", "V segments",
                    "J alleles", "J segments",
                    "D alleles", "D segments",
                    "Last V nucleotide position", "First D nucleotide position",
                    "Last D nucleotide position", "First J nucleotide position",
                    "VD insertions", "DJ insertions", "Total insertions"].join("\t"))

            clusterizedCloneSet.clones.each { Clone clone ->
                def nnnMap = new HashMap<String, Integer>()
                clone.backwardLinks.each { SequencingReadLink link ->
                    def nnnSeq = nnnList[(int) link.id]
                    nnnMap.put(nnnSeq, (nnnMap.get(nnnSeq) ?: 0) + 1)
                }
                int nnnTotal = 0
                nnnMap.values().each { if (it >= nnnCountThreshold) nnnTotal++ }
                pw.println([nnnTotal, clone.count, clone.count / clusterizedCloneSet.totalCount,
                        clone.CDR3.sequence, clone.CDR3.quality,
                        clone.CDR3AA,
                        writeAlleles(clone.VAlleles), writeSegments(clone.VSegments),
                        writeAlleles(clone.JAlleles), writeSegments(clone.JSegments),
                        writeAlleles(clone.DAlleles), writeSegments(clone.DSegments),
                        clone.VEnd, clone.DStart, clone.DEnd, clone.JStart,
                        clone.insertionsVD(), clone.insertionsDJ(), clone.insertionsTotal()
                ].join("\t"))
            }
        }
    }
})
new Thread(new ProgressReporter(pipeline)).start()
pipeline.run()