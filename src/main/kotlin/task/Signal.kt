package task

import krews.core.WorkflowBuilder
import krews.file.File
import krews.file.OutputFile
import model.*
import org.reactivestreams.Publisher

data class SignalParams(
    val smoothingWidth: Double = 30.0
)

data class SignalInput(
    val repName: String,
    val bams: List<File>
)

data class SignalOutput(
    val repName: String,
    val bothStrands: File,
    val plusStrand: File,
    val minusStrand: File
)

fun WorkflowBuilder.signalTask(i: Publisher<SignalInput>) = this.task<SignalInput, SignalOutput>("signal", i) {
    val params = taskParams<SignalParams>()

    dockerImage = "genomealmanac/mnase-signal:1.0.1"

    output = SignalOutput(
        repName = input.repName,
	bothStrands = OutputFile("signal/${input.repName}.wig"),
	plusStrand = OutputFile("signal/${input.repName}.plus.wig"),
	minusStrand = OutputFile("signal/${input.repName}.minus.wig")
    )
    
    command =
            """
            mkdir -p $dockerDataDir/signal && \
            mnasesignal . ${params.smoothingWidth} $dockerDataDir/signal/${input.repName}.wig ${input.bams.map{ it.dockerPath }.joinToString(" ")} && \
            mnasesignal + ${params.smoothingWidth} $dockerDataDir/signal/${input.repName}.plus.wig ${input.bams.map{ it.dockerPath }.joinToString(" ")} && \
            mnasesignal - ${params.smoothingWidth} $dockerDataDir/signal/${input.repName}.minus.wig ${input.bams.map{ it.dockerPath }.joinToString(" ")}
            """
}
