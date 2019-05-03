import krews.core.*
import krews.run
import model.*
import reactor.core.publisher.*
import task.*

fun main(args: Array<String>) = run(mnaseWorkflow, args)

data class MNaseParams(
    val samples: Samples
)

val mnaseWorkflow = workflow("mnase-workflow") {

    val params = params<MNaseParams>()

    val trimAdaptorInputs =
	if (params.samples is FastqSamples) {
	    params.samples.replicates
		.map { TrimAdaptorInput(it) }
		.toFlux()
        } else {
	    val fastaToFastqInputs = params.samples.replicates
	        .map { FastaToFastqInput(it as FastaReplicate) }
		.toFlux()
	    val fastaToFastqOutputs = fastaToFastqTask(fastaToFastqInputs)
	    fastaToFastqOutputs.map { TrimAdaptorInput(it.rep) }
        }
    val trimAdaptorOutputs = trimAdaptorTask(trimAdaptorInputs)

    val bowtie2Input = trimAdaptorOutputs
            .map { Bowtie2Input(it.mergedReplicate) }
    val bowtie2Task = bowtie2Task(bowtie2Input)

}
