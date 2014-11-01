package mitcr

@Grab(group = 'com.milaboratory', module = 'mitcr-groovy', version = '1.0.3')
@Grab(group = 'com.milaboratory', module = 'micommons', version = '1.0.3')
@Grab(group = 'cc.redberry', module = 'pipe', version = '0.9.2')

import cc.redberry.pipe.util.CountLimitingOutputPort
import com.milaboratory.core.clone.Clone
import com.milaboratory.core.clone.CloneSetClustered
import com.milaboratory.core.io.CloneSetIO
import com.milaboratory.core.segment.AlleleSet
import com.milaboratory.core.segment.SegmentSet
import com.milaboratory.core.sequencing.io.fastq.SFastqReader
import com.milaboratory.core.sequencing.read.SSequencingRead
import com.milaboratory.mitcr.cli.ExportDetalizationLevel
import com.milaboratory.mitcr.clonegenerator.SequencingReadLink
import com.milaboratory.mitcr.pipeline.DefaultAnalysisListener
import com.milaboratory.mitcr.pipeline.FullPipeline

/**
 INTERNAL
 */
import com.milaboratory.mitcr.pipeline.ParameterPresets
import com.milaboratory.util.CompressionType
import com.milaboratory.util.ProgressReporter

def params = ParameterPresets.jPrimer
def nReads = 10000
def reads = new SFastqReader("test_cdr.fastq.gz", CompressionType.GZIP)
reads = new CountLimitingOutputPort<SSequencingRead>(reads, nReads)
def pipeline = new FullPipeline(reads, params, true)
def header = ["CDR3 nucleotide sequence",
        "V alleles", "V segments",
        "J alleles", "J segments",
        "D alleles", "D segments",
        "Last V nucleotide position", "First D nucleotide position", "Last D nucleotide position", "First J nucleotide position",
        "VD insertions", "DJ insertions", "Total insertions"]
def writeAlleles = { AlleleSet alleles ->
    alleles.collect { it.fullName }.join(",")
}
def writeSegments = { SegmentSet segments ->
    segments.collect { it.segmentName }.join(",")
}
pipeline.setAnalysisListener(new DefaultAnalysisListener() {
    @Override
    void afterClusterization(CloneSetClustered clusterizedCloneSet) {
        def readMap = new HashMap<Integer, String>()
        for (Clone clone : clusterizedCloneSet.clones)
            for (SequencingReadLink readLink : clone.backwardLinks)
                readMap.put((int) readLink.id, [clone.CDR3.sequence,
                        writeAlleles(clone.VAlleles), writeSegments(clone.VSegments),
                        writeAlleles(clone.JAlleles), writeSegments(clone.JSegments),
                        writeAlleles(clone.DAlleles), writeSegments(clone.DSegments),
                        clone.VEnd, clone.DStart, clone.DEnd, clone.JStart,
                        clone.insertionsVD(), clone.insertionsDJ(), clone.insertionsTotal()
                ].join("\t"))
        new File("read_links.txt").withPrintWriter { pw ->
            pw.println("Read Id\t" + header.join("\t"))
            (0..(nReads - 1)).each {
                def value = readMap.get(it)
                pw.println(it + "\t" + value ?: header.collect { "?" }.join("\t"))
            }
        }
    }
})
new Thread(new ProgressReporter(pipeline)).start();
pipeline.run()

CloneSetIO.exportCloneSet("result.txt", pipeline.result, ExportDetalizationLevel.Full, params,
        "test_cdr.fastq.gz", CompressionType.None);