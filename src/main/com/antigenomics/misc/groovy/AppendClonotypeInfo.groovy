@Grab(group = 'com.milaboratory', module = 'mitcr-groovy', version = '1.0.3')
@Grab(group = 'com.milaboratory', module = 'micommons', version = '1.0.3')
@Grab(group = 'cc.redberry', module = 'pipe', version = '0.9.2')

import com.milaboratory.core.clone.Clone
import com.milaboratory.core.clone.CloneSetClustered
import com.milaboratory.core.segment.Gene
import com.milaboratory.core.segment.Species
import com.milaboratory.core.sequence.NucleotideSQPair
import com.milaboratory.core.sequence.quality.QualityFormat
import com.milaboratory.core.sequencing.io.fastq.SFastqReader
import com.milaboratory.core.sequencing.io.fastq.SFastqWriter
import com.milaboratory.core.sequencing.read.SSequencingReadImpl
import com.milaboratory.mitcr.clonegenerator.SequencingReadLink
import com.milaboratory.mitcr.pipeline.DefaultAnalysisListener
import com.milaboratory.mitcr.pipeline.FullPipeline
import com.milaboratory.mitcr.pipeline.ParameterPresets
import com.milaboratory.mitcr.qualitystrategy.IlluminaQualityInterpretationStrategy
import com.milaboratory.util.CompressionType
import com.milaboratory.util.ProgressReporter

if (args.length < 2) {
    println "groovy AppendClonotypeInfo input_R2 output_prefix MODE SPECIES GENE QUALITY"
    println "Allowed values: MODE=1 - will output a table with read ids and clonotypes, " +
            "2 - will append clonotype info to FASTQ header, 3 in[t,f]ernal mode, " +
            "SPECIES=HomoSapines,MusMusculus GENE=TRA,TRB QUALITY=2..40"
    println "Default values: MODE=1 SPECIES=HomoSapines GENE=TRB QUALITY=20"
    System.exit(0)
}

def r2FileName = args[0], outFileName = args[1]
def mode = Integer.parseInt(args.length > 2 ? args[2] : "1")
def species = args.length > 3 ? args[3] : "HomoSapiens", gene = args.length > 4 ? args[4] : "TRB"
def quality = Integer.parseInt(args.length > 5 ? args[5] : '20')

switch (mode) {
    case 1:
        outFileName += ".txt"
        break
    case 2:
    case 3:
        outFileName += r2FileName.endsWith(".fastq.gz") ? ".fastq.gz" : ".fastq"
        break
    default:
        println "Bad mode specified, $mode"
        System.exit(-1)
}

def params = ParameterPresets.jPrimer
params.species = Species."$species"
params.gene = Gene."$gene"
params.qualityInterpretationStrategy = new IlluminaQualityInterpretationStrategy((byte) quality)

def reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
def pipeline = new FullPipeline(reads, params, true)
def readMap = new HashMap<Long, String[]>()
pipeline.setAnalysisListener(new DefaultAnalysisListener() {
    @Override
    void afterClusterization(CloneSetClustered clusterizedCloneSet) {
        for (Clone clone : clusterizedCloneSet.clones)
            for (SequencingReadLink readLink : clone.backwardLinks) {
                def vIter = clone.VSegments.iterator(), jIter = clone.JSegments.iterator(),
                    dIter = clone.DSegments.iterator()
                readMap.put(readLink.id,
                        [
                                clone.CDR3.sequence,
                                vIter.hasNext() ? vIter.next().segmentName : "",
                                jIter.hasNext() ? jIter.next().segmentName : "",
                                dIter.hasNext() ? dIter.next().segmentName : ""
                        ] as String[])
            }

        if (mode == 1) {
            new File(outFileName).withPrintWriter { pw ->
                pw.println("#read_id\tcdr3\tv\tj")
                readMap.each {
                    pw.println([it.key, it.value.collect()].flatten().join("\t"))
                }
            }
        }
    }
})
new Thread(new ProgressReporter(pipeline)).start()
pipeline.run()

if (mode > 1) {
    reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
    def writer = new SFastqWriter(outFileName, QualityFormat.Phred33, outFileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)

    def read
    while ((read = reads.take()) != null) {
        def readInfo = readMap[read.id()]
        if (readInfo != null) {
            def descr = read.description + " CDR3:" + readInfo[0] + " V:" + readInfo[1] + " J:" + readInfo[2] + " D:" + readInfo[3]
            def data = mode == 2 ? read.getData() : new NucleotideSQPair("")
            def read2 = new SSequencingReadImpl(descr, data, read.id())
            writer.write(read2)
        }
    }
    writer.close()
}