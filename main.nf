#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/crisprvar
========================================================================================
    Github : https://github.com/nf-core/crisprvar
    Website: https://nf-co.re/crisprvar
    Slack  : https://nfcore.slack.com/channels/crisprvar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CRISPRVAR } from './workflows/crisprvar'

//
// WORKFLOW: Run main nf-core/crisprvar analysis pipeline
//
workflow NFCORE_CRISPRVAR {
    CRISPRVAR ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CRISPRVAR ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
