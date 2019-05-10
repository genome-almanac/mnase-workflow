package task

import krews.core.WorkflowBuilder
import krews.file.OutputFile
import model.*
import org.reactivestreams.Publisher

data class FastaToFastqInput(
    val rep: FastaReplicate
)

data class FastaToFastqOutput(
    val rep: FastqReplicate
)

fun WorkflowBuilder.fastaToFastqTask(i: Publisher<FastaToFastqInput>) = this.task<FastaToFastqInput, FastaToFastqOutput>("fasta-to-fastq", i) {

    dockerImage = "genomealmanac/mnase-fastatofastq:1.0.7"

    val rep = input.rep
    output = FastaToFastqOutput(
        FastqReplicateSE(
    	    name = rep.name,
	    adaptor = rep.adaptor,
	    fastqs = listOf(OutputFile("fastatofastq/${rep.name}.fastq"))
	)
    )

    command =
            """
            fastatofastq \
                ${rep.fasta.dockerPath} ${rep.csqual.dockerPath} $dockerDataDir/fastatofastq/${rep.name}.fastq
            """
}
