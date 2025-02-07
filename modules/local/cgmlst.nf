process chewbbacaAlleleCall {
    label "chewbbaca"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(assembly)
        path(schema_database)
    output:
        tuple val(meta), path("*.tsv")
    script:
    """
    chewBBACA.py AlleleCall -i ${assembly} -g ${schema_database} -o chewBBACA --no-inferred
    mv chewBBACA/results_alleles.tsv ${meta.alias}_results_alleles.tsv
    """
}

process chewbbacaJoinProfiles {
    label "chewbbaca"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(profile_1), path(profile_2)
    output:
        tuple val(meta), path("*.tsv")
    script:
    """
    chewBBACA.py JoinProfiles -p ${profile_1} ${profile_2} -o ${meta.alias}_joined_profile.tsv
    """
}

process grapetree {
    label "grapetree"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(joined_profile)
    output:
        tuple val(meta), path("*.nw")
    script:
    """
    grapetree -p ${joined_profile} > ${meta.alias}_cgMLST.nw
    """
}

workflow run_cgmlst {
   take:
        consensus
        scheme_database
        cgmlst_profile
   main:
        scheme = Channel.fromPath(scheme_database)
        cgmlstProfile = Channel.fromPath(cgmlst_profile)
        consensus_dir = consensus | map { meta, consensus -> [ meta, consensus.getParent() ] }
        // Alelle calling
        cgmlst_result = chewbbacaAlleleCall(consensus_dir, scheme.collect())
        profile_for_join = cgmlst_result | combine (cgmlstProfile)
        // cgMLST profile joining
        concat_cgmlst = chewbbacaJoinProfiles(profile_for_join)
        result_tree = grapetree(concat_cgmlst)
    emit:
        cgmlst_grapetree = result_tree
}