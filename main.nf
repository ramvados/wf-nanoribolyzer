#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                     Checkup imported Variables                                                          //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                     Basecalling Dorado                                                                  //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process dorado_basecalling{
    label 'dorado_basecaller'
    publishDir "${params.out_dir}/shared/", mode: 'copy'
    stageInMode 'symlink'
    input:
    path(sample_folder)
    val basecalling_model

    output:
    val 1, emit: done
    path('basecalling_output/basecalled.bam'), emit: basecalled_bam
    path('basecalling_output/sequencing_summary.txt'), optional: true
    path('basecalling_output/basecalled_not_trimmed.fastq.gz'), emit: fastq_not_trimmed
    path('converted_to_pod5/converted.pod5'), emit: converted_pod5
    
    script:
    def dorado_executable = params.get('dorado_executable') ?: 'dorado'
    def dorado_models_directory = params.get('dorado_models_directory') ?: ''
    """ 
    DORADO_MODELS_ARGS=()
    if [ -n "${dorado_models_directory}" ]; then
        DORADO_MODELS_ARGS=(--models-directory "${dorado_models_directory}")
    fi

    mkdir -p basecalling_output
    mkdir -p converted_to_pod5

    if compgen -G "${sample_folder}/*.pod5" > /dev/null; then
        export filetype=pod5
    elif compgen -G "${sample_folder}/*.fast5" > /dev/null; then
        export filetype=fast5
    else
        export filetype=bam
    fi 
    ################################
    if [ \$filetype == fast5 ]
    then 
        pod5 convert fast5 ${sample_folder}/*.fast5\
         --output converted_to_pod5/converted.pod5\
         --force-overwrite
        ${dorado_executable} basecaller "\${DORADO_MODELS_ARGS[@]}" ${basecalling_model} converted_to_pod5/\
         > basecalling_output/basecalled.bam 
        ${dorado_executable} summary basecalling_output/basecalled.bam\
         > basecalling_output/sequencing_summary.txt
        samtools bam2fq basecalling_output/basecalled.bam\
         -@ ${params.threads}\
         > basecalling_output/basecalled_not_trimmed.fastq
        gzip basecalling_output/basecalled_not_trimmed.fastq\
         -c\
         -1\
          > basecalling_output/basecalled_not_trimmed.fastq.gz
    fi
    if [ \$filetype == pod5 ]
    then
        ${dorado_executable} basecaller "\${DORADO_MODELS_ARGS[@]}" ${basecalling_model} ${sample_folder}\
         > basecalling_output/basecalled.bam
        ${dorado_executable} summary basecalling_output/basecalled.bam\
         > basecalling_output/sequencing_summary.txt
        samtools bam2fq basecalling_output/basecalled.bam\
         -@ ${params.threads}\
         > basecalling_output/basecalled_not_trimmed.fastq
        gzip basecalling_output/basecalled_not_trimmed.fastq\
         -c\
         -1\
          > basecalling_output/basecalled_not_trimmed.fastq.gz
        echo "No conversion needed" > converted_to_pod5/converted.pod5
    fi
    ###############################
    if [ \$filetype == bam ]
    then
        samtools view --threads ${params.threads} -bh ${sample_folder}/*bam > basecalling_output/basecalled.bam
        samtools fastq -T "*" --threads ${params.threads} basecalling_output/basecalled.bam | gzip > basecalling_output/basecalled_not_trimmed.fastq.gz
        echo "No conversion needed" > converted_to_pod5/converted.pod5
    fi
    """
}


process trim_barcodes{
    label 'barbell_tools'
    publishDir "${params.out_dir}/shared/barbell_demultiplexing/", mode: 'copy'
    stageInMode 'symlink'

    input:
        path(fastq_not_trimmed)
        path(barcode_fasta)
        path(barbell_filters)
        val(demultiplex)

    output:
        path("trimmed/*.fastq.gz"), emit: demultiplexed_fastqs, optional: true
        path("trimmed"), emit: trimmed_directory
        path("filtered.tsv"), emit: filter_table

    script:
    """
    zcat ${fastq_not_trimmed} > reads.fastq

    barbell annotate \
        -q ${barcode_fasta} \
        -b Ftag \
        -i reads.fastq \
        -o anno.tsv \
        -t ${params.threads}

    barbell filter \
        -i anno.tsv \
        -f ${barbell_filters} \
        -o filtered.tsv

    mkdir -p trimmed

    barbell trim \
        -i filtered.tsv \
        -r reads.fastq \
        -o trimmed \
        --gzip

    """
}


process retain_barbell_reads {
    label 'other_tools'
    publishDir "${params.out_dir}/shared/barbell_retention/", mode: 'copy'

    input:
        path(raw_fastq)
        path(trimmed_directory), stageAs: "barbell_trimmed"
        path(filter_table), stageAs: "filtered.tsv"
        path(fallback_script), stageAs: "barbell_fallback.py"
        path(context_barcode_fasta), stageAs: "context_barcodes.fasta"
        val(trim_fallback_context)

    output:
        path("basecalled.fastq.gz"), emit: combined_fastq
        path("samples"), emit: sample_directory
        path("retention_counts.tsv"), emit: counts

    script:
    """
    mkdir samples

    if compgen -G "barbell_trimmed/*.fastq.gz" > /dev/null; then
        cp barbell_trimmed/*.fastq.gz samples/
    fi

    python3 barbell_fallback.py \
        --raw ${raw_fastq} \
        --trimmed barbell_trimmed \
        --filtered filtered.tsv \
        --output samples/unassigned.fastq.gz \
        --counts retention_counts.tsv \
        --barcode-fasta context_barcodes.fasta \
        ${trim_fallback_context ? '--trim-context' : ''}

    if awk -F '\\t' '\$1 == "fallback_raw" && \$2 == 0 { found=1 } END { exit !found }' retention_counts.tsv; then
        rm samples/unassigned.fastq.gz
    fi

    zcat samples/*.fastq.gz | gzip -c > basecalled.fastq.gz
    """
}

process porechop_trimming{
    label 'other_tools'
    publishDir "${params.out_dir}/shared/porechop/", mode: 'copy'
    stageInMode 'symlink'

    input:
        path(fastq_not_trimmed)

    output:
        path("basecalled.fastq.gz"), emit: basecalled_fastq

    script:
    """
    mkdir -p chunks
    zcat ${fastq_not_trimmed} | split -l 8000000 -d -a 4 - chunks/chunk_
    : > basecalled.fastq

    for chunk in chunks/chunk_*; do
        porechop -i "\$chunk" -o "\${chunk}_trimmed" --threads ${params.threads}
        cat "\${chunk}_trimmed" >> basecalled.fastq
        rm "\$chunk" "\${chunk}_trimmed"
    done

    gzip basecalled.fastq
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                        Alignment of basecalled sequence to reference seq                                                //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process align_to_45SN1{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/basecalling_output/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(basecalled_fastq)
        path(reference), stageAs:"reference.fasta"
    output:
        tuple val(sample_id), path("filtered.bam"), path("filtered.bam.bai"), emit: aligned_bam
        tuple val(sample_id), path("filtered.fastq.gz"), emit: filtered_fastq
        tuple val(sample_id), val(1), emit: done
    script:
    """
    minimap2\
     -ax map-ont\
     -t ${params.threads}\
     --MD\
     reference.fasta\
     ${basecalled_fastq}\
     | samtools view -hbS -F 3884\
     | samtools sort\
     > filtered.bam
    #################
    samtools bam2fq filtered.bam --threads ${params.threads} > filtered.fastq
    ##################
    samtools index filtered.bam -@ ${params.threads}
    gzip filtered.fastq -c > filtered.fastq.gz
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                             Filter for reads in fast5 that align 45SN1                                                  //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process filter_pod5_for_RNA45s_aligning_reads{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/filtered_pod5/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(filtered_bam), path(filtered_bai)
        path(sample_folder)
        path(converted_pod5), stageAs: "converted.pod5"
    output:
        tuple val(sample_id), path("filtered.pod5"), path("filtered.bam"), emit: filtered_reads
        tuple val(sample_id), path("sorted_filtered_reads.txt"), emit: sorted_filtered_read_ids, optional: true
        tuple val(sample_id), val(1), emit: done
    script:
    """
    mkdir -p filtered_pod5
    (ls ${sample_folder}/*.pod5) && export filetype=pod5 || export filetype=fast5 
    
    if compgen -G "${sample_folder}/*.pod5" > /dev/null; then
        export filetype=pod5
    elif compgen -G "${sample_folder}/*.fast5" > /dev/null; then
        export filetype=fast5
    else
        export filetype=bam
    fi 
    echo \$filetype
    
    if [ \$filetype == pod5 ]
    then
        python ${projectDir}/bin/filter_pod5.py\
         -i ${filtered_bam}\
         -o .
        pod5 filter ${sample_folder}/*.pod5\
         --ids sorted_filtered_reads.txt\
         --output filtered.pod5\
         --missing-ok\
         --force-overwrite\
         --threads ${params.threads}
    fi
    if [ \$filetype == fast5 ]
    then
        python ${projectDir}/bin/filter_pod5.py\
         -i ${filtered_bam}\
         -o .
        pod5 filter converted.pod5\
         --ids sorted_filtered_reads.txt\
         --output filtered.pod5\
         --missing-ok\
         --force-overwrite\
         --threads ${params.threads}
    fi
    if [ \$filetype == bam ]
    then
        echo "" > filtered.pod5
    fi
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Rebasecalling in order to obtain fast5_out option and Event tables in fast5                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process rebasecall_filtered_files{
    label 'dorado_basecaller'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/filtered_pod5/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(filtered_pod5), path(filtered_bam)
        path(reference), stageAs: "reference.fasta"
        val basecalling_model
        val filetype
    output:
        tuple val(sample_id), path("filtered_pod5_basecalled.bam"), path("filtered_pod5_basecalled.bam.bai"), emit: rebasecalled
        tuple val(sample_id), val(1), emit: done
    script:
    def dorado_executable = params.get('dorado_executable') ?: 'dorado'
    def dorado_models_directory = params.get('dorado_models_directory') ?: ''
    def dorado_rna_model = params.get('dorado_rna_model') ?: 'sup,m6A,pseU'
    def dorado_modified_bases_models = params.get('dorado_modified_bases_models') ?: ''
    """
    DORADO_MODELS_ARGS=()
    if [ -n "${dorado_models_directory}" ]; then
        DORADO_MODELS_ARGS=(--models-directory "${dorado_models_directory}")
    fi

    DORADO_MOD_ARGS=()
    if [ -n "${dorado_modified_bases_models}" ]; then
        DORADO_MOD_ARGS=(--modified-bases-models "${dorado_modified_bases_models}")
    fi

    if [ ${filetype} == "bam" ]
    then
        mv ${filtered_bam} filtered_pod5_basecalled.bam
        samtools index filtered_pod5_basecalled.bam
    else
        if [ ${params.sample_type} == "DNA" ]
        then
            ${dorado_executable} basecaller "\${DORADO_MODELS_ARGS[@]}" --estimate-poly-a --emit-moves ${basecalling_model} ${filtered_pod5}\
            | samtools fastq -T "*" --threads ${params.threads}\
            | minimap2 -t ${params.threads} -y --MD -ax map-ont reference.fasta -\
            | samtools sort --threads ${params.threads}\
            | samtools view -b -F 3884 --threads ${params.threads}\
            > filtered_pod5_basecalled.bam 
            samtools index filtered_pod5_basecalled.bam -@ ${params.threads}
        else
            ${dorado_executable} basecaller --device "cuda:0" --estimate-poly-a --emit-moves "\${DORADO_MOD_ARGS[@]}" ${dorado_rna_model} ${filtered_pod5}\
            | samtools fastq -T "*" --threads ${params.threads}\
            | minimap2 -t ${params.threads} -y --MD -ax map-ont reference.fasta -\
            | samtools sort --threads ${params.threads}\
            | samtools view -b -F 3884 --threads ${params.threads}\
            > filtered_pod5_basecalled.bam 
            samtools index filtered_pod5_basecalled.bam -@ ${params.threads}
        fi
    fi
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                          Run polyA estimation based on tailfindr                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process extract_polyA_table{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/taillength_estimation/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(rebasecalled_bam), path(rebasecalled_bam_bai)
    output:
        tuple val(sample_id), path("tail_estimation.csv"), emit: tail_estimation_csv
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/extract_polyA_tails.py -i filtered_pod5_basecalled.bam -o .
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis without template to obtain sample specific fragment cluster                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis without template to obtain sample specific fragment cluster                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process fragment_analysis_intensity{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/fragment_analysis_intensity/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(filtered_bam), path(filtered_bam_bai)
        path(reference), stageAs: "reference.fasta" 
        // val(done)
    output:
        tuple val(sample_id), path("alignment_df.csv"), emit: alignment_df
        tuple val(sample_id), path("fragment_df.csv"), emit: fragment_df
        tuple val(sample_id), path("fragment_analysis_intensity_matrix.png"), emit: intensity_matrix_png
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/fragment_analysis_intensity.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -o ./ -t 0.9 -s ${params.color} -d ${params.demand}
    python ${projectDir}/bin/visualize_clustering_performance_intensity_matrix.py -a ./alignment_df.csv -t ./fragment_df_simple.csv -c ${params.color} -o ./
    cp intensity_matrix.png fragment_analysis_intensity_matrix.png
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                           Run fragment analysis with templates known to be present in literature                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process template_driven_fragment_analysis{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/template_based_analysis/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(filtered_bam), path(filtered_bam_bai), val(done)
        path(reference), stageAs: "reference.fasta" 
        path(templates), stageAs: "template.csv"
    output:
        tuple val(sample_id), path("template_fragment_df.csv"), emit: template_csv
        tuple val(sample_id), path("template_alignment_df.csv"), emit: template_alignment_csv
        tuple val(sample_id), path("start_sites_fragment_based.bed"), emit: start_sites
        tuple val(sample_id), path("end_sites_fragment_based.bed"), emit: end_sites
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/template_driven_fragment_analysis.py -c ${params.threads} -i ${filtered_bam} -r reference.fasta -f template.csv -o ./ -t 0.9 -s ${params.color} -d ${params.demand}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                               Run read-tail analysis based on associated templates of template driven fragment analysis                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process fragment_based_readtail_analysis{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/readtail_analysis/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(filtered_bam), path(filtered_bam_bai), path(template_csv)
        path(fasta_file)
    output:
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/readtail_analysis.py -i ${filtered_bam} -r ${template_csv} -o ./ -f ${fasta_file}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for template base associations                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_polyA_associated_templates{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/polyA_template_based/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(template_alignment_csv), path(tail_estimation_csv)
        path(templates)
        path(reference)
    output:
        tuple val(sample_id), path("polyA_tails_intermediates_template.html"), path("polyA_tails_intermediates_min_max.html"), path("polyA_tails_intermediates_mean.html"), path("violinplot_taillength_per_intermediate.png"), emit: report_files
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/visualize_taillengths.py -a ${template_alignment_csv} -t ${tail_estimation_csv} -c ${params.color} -o ./ -f ${templates} -r ${reference} -m ${params.model_organism}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for instensity based clusters                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_polyA_associated_intensity_clusters{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/polyA_intensity_based_clusters/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(intensity_alignment_csv), path(tail_estimation_csv)
        path(templates)
        path(reference)
    output:
        tuple val(sample_id), path("polyA_tails_clustering.html"), path("violinplot_taillength_per_cluster.png"), emit: report_files
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/visualize_taillengths_clustering.py -a ${intensity_alignment_csv} -t ${tail_estimation_csv} -c ${params.color} -o ./ -f ${templates} -r ${reference} -m ${params.model_organism}
    """
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                      Visualize poly-A taillengths for instensity based clusters                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize and save intensity matrix                                                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_intensity_matrix{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/intensity_matrix/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(alignment_df)
        path(templates)
    output:
        tuple val(sample_id), path("intensity_matrix.png"), path("intensity_matrix.html"), emit: report_files
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/visualize_intensity_matrix.py -a ${alignment_df} -t ${templates} -c ${params.color} -o ./ -m ${params.model_organism}
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize modifications (PsU,m6A)                                                               //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_modifications{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/modification_plots/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(rebasecalled_bam), path(rebasecalled_bam_bai)
        path(reference)
        path(modifications_bed)
    output:
        tuple val(sample_id), path("relative_pseU_modification_abundance.html"), path("relative_m6A_modification_abundance.html"), emit: report_files, optional: true
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    if [ ${params.sample_type} == "RNA" ]
    then
        python ${projectDir}/bin/visualize_modifications.py -b ${rebasecalled_bam} -r ${reference} -m ${modifications_bed} -o ./
    else
        echo "cDNA used" > modifications.log
    fi
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize significant cut sites                                                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_cut_sites{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/cut_site_plots/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(start_sites), path(end_sites)
        path(reference)
    output:
        tuple val(sample_id), path("cut_sites.html"), emit: report_files
        tuple val(sample_id), path("*"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/visualize_cut_sites.py -s ${start_sites} -e ${end_sites}  -r ${reference} -o ./
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                         Visualize reference coverage                                                                    //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualize_reference_coverage{
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/coverage_plots/" }, mode:"copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), path(template_alignment_csv)
        path(templates)
        path(reference)
    output:
        tuple val(sample_id), path("coverage_fragments_absolute.png"), path("coverage_fragments_relative.png"), path("coverage_fragments_absolute_all.png"), path("coverage_total_sample_absolute.png"), path("coverage_total_sample_relative.png"), emit: report_files
        tuple val(sample_id), path("*.png"), emit: all_files
        tuple val(sample_id), val(1), emit: done
    script:
    """
    python ${projectDir}/bin/visualize_rRNA_coverage.py -a ${template_alignment_csv} -f ${templates}  -r ${reference} -o ./ -c ${params.color} -m ${params.model_organism}
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                               Check if all workflows are done                                                           //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process check_all_done {
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/report/" }, mode: "copy"
    stageInMode 'symlink'
    input:
        tuple val(sample_id), val(visualize_intensity_matrix_done), val(visualize_polyA_associated_templates_done), val(visualize_polyA_associated_intensity_clusters_done), val(visualize_modifications_done), val(visualize_cut_sites_done), val(visualize_reference_coverage_done)
    output:
        tuple val(sample_id), val(1), emit: done
    script:
    """
    echo "All visualizations done"
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                               Create HTML report                                                                        //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process create_report {
    label 'other_tools'
    tag "$sample_id"
    publishDir { "${params.out_dir}/${sample_id}/" }, mode: "copy"
    stageInMode 'symlink'
    input:
    // Mandatory inputs
    tuple val(sample_id), val(all_done_verification), path(intensity_matrix_png), path(intensity_matrix_html), path(polyA_tails_intermediates_template_html), path(polyA_tails_intermediates_min_max_html), path(polyA_tails_intermediates_mean_html), path(violinplot_taillength_per_intermediate_png), path(cut_sites_html), path(polyA_tails_clustering_html), path(intensity_matrix_intensity_clustering), path(violinplot_taillength_per_cluster_png), path(coverage_plot_general_absolute), path(coverage_plot_general_relative), path(coverage_plot_fragments_absolute), path(coverage_total_sample_absolute), path(coverage_total_sample_relative), path(relative_pseU_modification_abundance_html), path(relative_m6A_modification_abundance_html)
    output:
    tuple val(sample_id), path("rRNA_report.html"), emit: report
    
    script:
    """
    mkdir -p intensity_matrix polyA_template_based cut_site_plots polyA_intensity_based_clusters fragment_analysis_intensity coverage_plots modification_plots
    cp ${intensity_matrix_png} intensity_matrix/intensity_matrix.png
    cp ${intensity_matrix_html} intensity_matrix/intensity_matrix.html
    cp ${polyA_tails_intermediates_template_html} polyA_template_based/polyA_tails_intermediates_template.html
    cp ${polyA_tails_intermediates_min_max_html} polyA_template_based/polyA_tails_intermediates_min_max.html
    cp ${polyA_tails_intermediates_mean_html} polyA_template_based/polyA_tails_intermediates_mean.html
    cp ${violinplot_taillength_per_intermediate_png} polyA_template_based/violinplot_taillength_per_intermediate.png
    cp ${cut_sites_html} cut_site_plots/cut_sites.html
    cp ${polyA_tails_clustering_html} polyA_intensity_based_clusters/polyA_tails_clustering.html
    cp ${intensity_matrix_intensity_clustering} fragment_analysis_intensity/intensity_matrix.png
    cp ${violinplot_taillength_per_cluster_png} polyA_intensity_based_clusters/violinplot_taillength_per_cluster.png
    cp ${coverage_plot_general_absolute} coverage_plots/coverage_fragments_absolute.png
    cp ${coverage_plot_general_relative} coverage_plots/coverage_fragments_relative.png
    cp ${coverage_plot_fragments_absolute} coverage_plots/coverage_fragments_absolute_all.png
    cp ${coverage_total_sample_absolute} coverage_plots/coverage_total_sample_absolute.png
    cp ${coverage_total_sample_relative} coverage_plots/coverage_total_sample_relative.png
    cp ${relative_pseU_modification_abundance_html} modification_plots/relative_pseU_modification_abundance.html
    cp ${relative_m6A_modification_abundance_html} modification_plots/relative_m6A_modification_abundance.html
    python ${projectDir}/bin/html_report.py -d ./ -o ./
    """
}



workflow{
    println "Seqfolder"
    println params.sample_folder
    println ""
    println "Color"
    println params.color
    println ""
    println "Script folder"
    println params.script_folder
    println ""
    println "Output folder"
    println params.out_dir
    println ""
    println "Basecalling model"
    println params.basecalling_model
    println ""
    println "Threads"
    println params.threads
    println ""
    println "Sample type"
    println params.sample_type
    println ""
    println "Demand"
    println params.demand
    println ""
    println "Model Organism"
    println params.model_organism
    println ""
    println "Demultiplex samples"
    println(params.demultiplex ?: false)
    println ""
    println "Selected barcodes"
    println(params.get('barcodes') ?: "all detected barcodes")
    println ""

    if (params.model_organism == "Human"){
        fasta_reference_file = "${projectDir}/references/RNA45SN1.fasta"
        ribosomal_intermediates_file = "${projectDir}/references/Literature_Fragments_and_cut_sites_RNA45SN1.csv"
        modification_reference_file = "${projectDir}/references/rRNA_modifications_conv.bed"
    } else if (params.model_organism == "Yeast"){
        fasta_reference_file = "${projectDir}/references/RDN37-1.fa"
        ribosomal_intermediates_file = "${projectDir}/references/Literature_Fragments_and_cut_sites_RDN37-1.csv"
        modification_reference_file = "${projectDir}/references/rRNA_yeast_modifications_conv.bed"
    } else {
        error "Unsupported model_organism: ${params.model_organism}. Use Human or Yeast."
    }

    println fasta_reference_file
    println ribosomal_intermediates_file
    println modification_reference_file

    def sample_dir = file(params.sample_folder)
    def filetype = sample_dir.listFiles().any { entry -> entry.name.endsWith('.pod5') } ? 'pod5' :
                   sample_dir.listFiles().any { entry -> entry.name.endsWith('.fast5') } ? 'fast5' :
                   sample_dir.listFiles().any { entry -> entry.name.endsWith('.bam') } ? 'bam' : null
    if (!filetype) {
        error "No POD5, FAST5 or BAM input found in ${params.sample_folder}"
    }

    def preprocessing_method = (params.preprocessing_method ?: 'barbell').toString().trim().toLowerCase()
    if (!['porechop', 'barbell'].contains(preprocessing_method)) {
        error "Unsupported preprocessing_method: ${preprocessing_method}. Use porechop or barbell."
    }

    def demultiplex_enabled = params.demultiplex != null && params.demultiplex.toString().toBoolean()
    if (preprocessing_method == 'porechop' && demultiplex_enabled) {
        error "demultiplex=true requires preprocessing_method=barbell"
    }

    def selected_barcodes = params.get('barcodes') ?: []
    if (selected_barcodes instanceof String) {
        selected_barcodes = selected_barcodes.split(',')*.trim().findAll { barcode -> barcode }
    }

    dorado_basecalling(
        file(params.sample_folder),
        params.basecalling_model
    )


    if (preprocessing_method == 'barbell') {
        trim_barcodes(
            dorado_basecalling.out.fastq_not_trimmed,
            file(params.containsKey('barbell_barcode_fasta') && params.barbell_barcode_fasta ? params.barbell_barcode_fasta : "${projectDir}/data/DRB004_RNA01-12.fasta"),
            file(params.containsKey('barbell_filters') && params.barbell_filters ? params.barbell_filters : "${projectDir}/data/barbell_DRB004_filters.txt"),
            demultiplex_enabled
        )

        retain_barbell_reads(
            dorado_basecalling.out.fastq_not_trimmed,
            trim_barcodes.out.trimmed_directory,
            trim_barcodes.out.filter_table,
            file("${projectDir}/bin/barbell_fallback.py"),
            file(params.containsKey('barbell_barcode_fasta') && params.barbell_barcode_fasta ? params.barbell_barcode_fasta : "${projectDir}/data/DRB004_RNA01-12.fasta"),
            params.containsKey('barbell_trim_fallback_context') && params.barbell_trim_fallback_context != null
                ? params.barbell_trim_fallback_context.toString().toBoolean()
                : false
        )

        if (demultiplex_enabled) {
            sample_fastqs = retain_barbell_reads.out.sample_directory
                .flatMap { directory ->
                    directory.toFile().listFiles()
                        .findAll { entry -> entry.name.endsWith('.fastq.gz') }
                        .collect { entry -> file(entry.toPath()) }
                }
                .map { fastq ->
                    def sample_id = fastq.name
                        .replaceFirst(/\.fastq\.gz$/, '')
                        .replaceFirst(/_(fw|rc)\.trimmed$/, '')
                        .replaceFirst(/\.trimmed$/, '')
                    tuple(sample_id, fastq)
                }

            if (selected_barcodes) {
                sample_fastqs = sample_fastqs.filter { sample_id, _fastq ->
                    sample_id == 'unassigned' || selected_barcodes.contains(sample_id)
                }
            }
        } else {
            def sample_name = params.sample_name ?: 'sample'
            sample_fastqs = retain_barbell_reads.out.combined_fastq.map { fastq ->
                tuple(sample_name, fastq)
            }
        }
    } else {
        porechop_trimming(
            dorado_basecalling.out.fastq_not_trimmed
        )

        def sample_name = params.sample_name ?: 'sample'
        sample_fastqs = porechop_trimming.out.basecalled_fastq.map { fastq ->
            tuple(sample_name, fastq)
        }
    }

    sample_fastqs = sample_fastqs
        .ifEmpty { error "No FASTQ files were produced using ${preprocessing_method}; requested samples: ${selected_barcodes ?: 'all'}" }
        .view { sample_id, fastq -> "NanoRibolyzer sample [${preprocessing_method}]: ${sample_id} (${fastq.name})" }

    align_to_45SN1(sample_fastqs, file(fasta_reference_file))

    filter_pod5_for_RNA45s_aligning_reads(
        align_to_45SN1.out.aligned_bam,
        file(params.sample_folder),
        dorado_basecalling.out.converted_pod5
    )

    rebasecall_filtered_files(
        filter_pod5_for_RNA45s_aligning_reads.out.filtered_reads,
        file(fasta_reference_file),
        params.basecalling_model,
        filetype
    )

    extract_polyA_table(rebasecall_filtered_files.out.rebasecalled)
    fragment_analysis_intensity(
        rebasecall_filtered_files.out.rebasecalled,
        file(fasta_reference_file)
    )

    template_analysis_input = rebasecall_filtered_files.out.rebasecalled
        .join(fragment_analysis_intensity.out.done)
    template_driven_fragment_analysis(
        template_analysis_input,
        file(fasta_reference_file),
        file(ribosomal_intermediates_file)
    )

    readtail_input = rebasecall_filtered_files.out.rebasecalled
        .join(template_driven_fragment_analysis.out.template_csv)
    fragment_based_readtail_analysis(readtail_input, file(fasta_reference_file))

    visualize_intensity_matrix(
        template_driven_fragment_analysis.out.template_alignment_csv,
        file(ribosomal_intermediates_file)
    )

    polyA_template_input = template_driven_fragment_analysis.out.template_alignment_csv
        .join(extract_polyA_table.out.tail_estimation_csv)
    visualize_polyA_associated_templates(
        polyA_template_input,
        file(ribosomal_intermediates_file),
        file(fasta_reference_file)
    )

    polyA_intensity_input = fragment_analysis_intensity.out.fragment_df
        .join(extract_polyA_table.out.tail_estimation_csv)
    visualize_polyA_associated_intensity_clusters(
        polyA_intensity_input,
        file(ribosomal_intermediates_file),
        file(fasta_reference_file)
    )

    visualize_modifications(
        rebasecall_filtered_files.out.rebasecalled,
        file(fasta_reference_file),
        file(modification_reference_file)
    )

    cut_sites_input = template_driven_fragment_analysis.out.start_sites
        .join(template_driven_fragment_analysis.out.end_sites)
    visualize_cut_sites(cut_sites_input, file(fasta_reference_file))

    visualize_reference_coverage(
        template_driven_fragment_analysis.out.template_alignment_csv,
        file(ribosomal_intermediates_file),
        file(fasta_reference_file)
    )

    all_visualizations_done = visualize_intensity_matrix.out.done
        .join(visualize_polyA_associated_templates.out.done)
        .join(visualize_polyA_associated_intensity_clusters.out.done)
        .join(visualize_modifications.out.done)
        .join(visualize_cut_sites.out.done)
        .join(visualize_reference_coverage.out.done)
    check_all_done(all_visualizations_done)

    report_inputs = check_all_done.out.done
        .join(visualize_intensity_matrix.out.report_files)
        .join(visualize_polyA_associated_templates.out.report_files)
        .join(visualize_cut_sites.out.report_files)
        .join(visualize_polyA_associated_intensity_clusters.out.report_files)
        .join(fragment_analysis_intensity.out.intensity_matrix_png)
        .join(visualize_reference_coverage.out.report_files)
        .join(visualize_modifications.out.report_files)
    create_report(report_inputs)
}
