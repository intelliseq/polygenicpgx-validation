#!/usr/bin/env nextflow
/*
 * PGx tool benchmark — polygenic vs. the field on GeT-RM x 1000 Genomes samples.
 *
 *   nextflow run benchmark/main.nf -profile docker
 *   nextflow run benchmark/main.nf -profile docker --samples NA12878 --genes cyp2c19 --tools polygenic
 *
 * Stages: FETCH (gene-slice VCF/CRAM per build) -> RUN_TOOL (each caller) ->
 *         HARMONIZE (native output -> canonical calls.tsv) -> SCORE (vs truth) -> REPORT.
 */
nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

// ---- registries -------------------------------------------------------------
def genes = new Yaml().load(file(params.genes_yml).text).genes
def tools = new Yaml().load(file("${projectDir}/conf/tool_versions.yml").text).tools

def wanted(sel, key) { sel == 'all' || sel.toString().split(',').contains(key) }

// tool applies to (gene, build)?  Group C tools never run.
def toolApplies = { String tool, String gene, String build ->
    def t = tools[tool]
    if (!t || t.group == 'C') return false
    if (!(build in (t.builds ?: []))) return false
    if (t.genes && !(gene in t.genes)) return false              // gene-restricted (e.g. cyrius=cyp2d6)
    if (tool == 'polygenic' && genes[gene].polygenic_model == null) return false  // unmodelled -> N/A
    return true
}

workflow {
    // sample, population, longread
    samples = Channel.fromPath(params.samples_tsv)
        .splitCsv(header: true, sep: '\t')
        .filter { wanted(params.samples, it.sample) }

    geneList  = genes.keySet().findAll { wanted(params.genes, it) }
    buildList = params.builds.split(',') as List
    toolList  = tools.keySet().findAll { wanted(params.tools, it) && tools[it].group in ['A', 'B'] }

    // (sample, gene, build) grid -> FETCH inputs
    grid = samples.combine(Channel.fromList(geneList)).combine(Channel.fromList(buildList))
    fetched = FETCH(grid.map { s, g, b -> tuple(s.sample, g, b, s.longread) })

    // expand each fetched bundle across the applicable tools
    jobs = fetched.flatMap { sample, gene, build, longread, bundle ->
        toolList.findAll { toolApplies(it, gene, build) ||
                           (tools[it].group == 'B' && longread == 'yes' && (tools[it].genes == null || gene in tools[it].genes)) }
                .collect { tuple(it, tools[it].image, sample, gene, build, bundle) }
    }

    raw = RUN_TOOL(jobs)
    calls = HARMONIZE(raw.collect())
    conc  = SCORE(calls, file(params.truth_tsv))
    REPORT(conc, calls)
}

process FETCH {
    tag "${sample}/${gene}/${build}"
    maxForks 6                       // cap concurrent remote-tabix pulls (1000G FTP throttles)
    errorStrategy 'retry'
    maxRetries 2
    input:  tuple val(sample), val(gene), val(build), val(longread)
    output: tuple val(sample), val(gene), val(build), val(longread), path("bundle")
    script:
    """
    python3 ${projectDir}/bin/fetch_data.py \
        --sample ${sample} --gene ${gene} --build ${build} \
        --genes-yml ${params.genes_yml} \
        --longread ${longread} \
        --out bundle
    """
}

process RUN_TOOL {
    tag "${tool}:${sample}/${gene}/${build}"
    container "${image}"
    publishDir "${params.outdir}/raw", mode: 'copy'   // persist raws for re-scoring / checkpoints
    input:  tuple val(tool), val(image), val(sample), val(gene), val(build), path(bundle)
    output: path("${tool}-${sample}-${gene}-${build}.raw")
    script:
    // bin/ scripts are auto-staged onto PATH by Nextflow (also inside containers).
    """
    run_${tool}.sh \
        ${sample} ${gene} ${build} ${bundle} \
        ${tool}-${sample}-${gene}-${build}.raw \
        || printf 'tool\\t%s\\nsample\\t%s\\ngene\\t%s\\nbuild\\t%s\\nstatus\\tERROR\\n' \
            ${tool} ${sample} ${gene} ${build} > ${tool}-${sample}-${gene}-${build}.raw
    """
}

process HARMONIZE {
    input:  path(raw_files)
    output: path("calls.tsv")
    publishDir params.outdir, mode: 'copy'
    script:
    """
    python3 ${projectDir}/bin/harmonize.py --raw ${raw_files} --out calls.tsv
    """
}

process SCORE {
    input:  path(calls); path(truth)
    output: path("concordance.tsv")
    publishDir params.outdir, mode: 'copy'
    script:
    """
    python3 ${projectDir}/bin/score.py --calls ${calls} --truth ${truth} \
        --samples ${params.samples_tsv} --out concordance.tsv
    """
}

process REPORT {
    input:  path(concordance); path(calls)
    output: path("report.md"); path("capability_matrix.md")
    publishDir params.outdir, mode: 'copy'
    script:
    """
    python3 ${projectDir}/bin/report.py \
        --concordance ${concordance} --calls ${calls} \
        --tool-versions ${projectDir}/conf/tool_versions.yml \
        --genes-yml ${params.genes_yml} \
        --report report.md --matrix capability_matrix.md
    """
}
