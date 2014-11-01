/**
 Copyright 2013 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package mitcr

@Grab(group = 'com.milaboratory', module = 'mitcr-groovy', version = '1.0.3')
import com.milaboratory.core.clone.Clone
import com.milaboratory.core.clone.CloneSetClustered
import com.milaboratory.core.segment.Gene
import com.milaboratory.core.segment.Species
import com.milaboratory.core.sequence.quality.QualityFormat
import com.milaboratory.core.sequencing.io.fastq.SFastqReader
import com.milaboratory.core.sequencing.io.fastq.SFastqWriter
import com.milaboratory.mitcr.clonegenerator.BasicCloneGeneratorParameters
import com.milaboratory.mitcr.clonegenerator.LQFilteringOffCloneGeneratorParameters
import com.milaboratory.mitcr.clonegenerator.SequencingReadLink
import com.milaboratory.mitcr.clusterization.CloneClusterizationType
import com.milaboratory.mitcr.pipeline.DefaultAnalysisListener
import com.milaboratory.mitcr.pipeline.FullPipeline
import com.milaboratory.mitcr.pipeline.ParameterPresets
import com.milaboratory.mitcr.qualitystrategy.IlluminaQualityInterpretationStrategy
import com.milaboratory.util.CompressionType
import com.milaboratory.util.ProgressReporter

if (args.length < 2) {
    println "groovy ExtractUniqueGoodSeqs input_R2 output SPECIES GENE QUALITY"
    println "Extracts good sequences, ie. those that 1) have a CDR3 region 2) CDR3 min quality is above threshold"
    println "Samples only reads that have unique UMIs"
    println ""
    println "Allowed values: SPECIES=HomoSapines,MusMusculus GENE=TRA,TRB QUALITY=2..40"
    println "Default values: SPECIES=HomoSapines GENE=TRB QUALITY=20"
    System.exit(0)
}

def r2FileName = args[0], outFileName = args[1]
def species = args.length > 2 ? args[2] : "HomoSapiens", gene = args.length > 3 ? args[3] : "TRB"
def quality = Integer.parseInt(args.length > 4 ? args[4] : '20')

println "[${new Date()}] Running ExtractUniqueGoodSeqs for $r2FileName"

def nnnList = new ArrayList<String>()
def nnnSet = new HashSet<String>()
def reads

// Read UMIs from R1
reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
def read
while ((read = reads.take())) {
    def umi = read.getDescription().split(" UMI:")[1].split(":")[0]
    nnnList.add(umi)
    nnnSet.add(umi)
}
println "[${new Date()}] ${nnnList.size()} reads loaded, ${nnnSet.size()} unique events"

def params = ParameterPresets.flex
params.species = Species."$species"
params.gene = Gene."$gene"
params.qualityInterpretationStrategy = new IlluminaQualityInterpretationStrategy((byte) quality)

// Disable assembly - essential:
params.clusterizationType = CloneClusterizationType.None
params.setCloneGeneratorParameters(new LQFilteringOffCloneGeneratorParameters(
        ((BasicCloneGeneratorParameters) params.getCloneGeneratorParameters()).getSegmentInformationAggregationFactor()))

reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
def pipeline = new FullPipeline(reads, params, true)

def usedNNNSet = new HashSet<String>()
def backwardLinks = new HashMap<Integer, String>()
int nGoodSeqs = 0
println "[${new Date()}] Running MiTCR"
pipeline.setAnalysisListener(new DefaultAnalysisListener() {
    @Override
    void afterClusterization(CloneSetClustered clusterizedCloneSet) {
        clusterizedCloneSet.clones.each { Clone clone ->
            clone.backwardLinks.each { SequencingReadLink link ->
                def nnn = nnnList[(int) link.id]
                if (!usedNNNSet.contains(nnn)) {
                    backwardLinks.put((Integer) link.id, nnn)
                    usedNNNSet.add(nnn)
                    nGoodSeqs++
                }
            }
        }
    }
})
new Thread(new ProgressReporter(pipeline)).start()
pipeline.run()
println "[${new Date()}] Extracted $nGoodSeqs good sequences"

// Read from R2 and compare IDs
println "[${new Date()}] Writing output"
def writer = new SFastqWriter(outFileName, QualityFormat.Phred33,
        outFileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
reads = new SFastqReader(r2FileName, r2FileName.endsWith(".gz") ? CompressionType.GZIP : CompressionType.None)
Integer n = 0
while (read = reads.take()) {
    if (backwardLinks.get(n))
        writer.write(read)
    n++
}
writer.close()