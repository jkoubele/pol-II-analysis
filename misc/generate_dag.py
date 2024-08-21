import graphviz
import pydot

if __name__ == "__main__":
    graph = graphviz.Digraph('DAG', filename='dag.gv')

    umi_color = 'lightgreen'
    default_color = 'darkblue'

    graph.node('FASTQ_and_UMI', color=umi_color)
    graph.node('FASTQ', color=default_color)
    graph.edge('FASTQ_and_UMI', 'FASTQ', color=umi_color, label=' Add UMI to read names')

    graph.node('QC_before_trimming', color=default_color)
    graph.edge('FASTQ', 'QC_before_trimming', label='FASTQC')

    graph.node('detected_adapters', color=default_color)
    graph.edge('FASTQ', 'detected_adapters', label='Atria - detect adapters')

    graph.node('detected_adapters_aggregated', color=default_color)
    graph.edge('detected_adapters', 'detected_adapters_aggregated', label='aggregate adapters')

    graph.node('FASTQ_trimmed', color=default_color)
    graph.edge('detected_adapters_aggregated', 'FASTQ_trimmed')
    graph.edge('FASTQ', 'FASTQ_trimmed', label='Atria trimming')

    graph.node('QC_after_trimming', color=default_color)
    graph.edge('FASTQ_trimmed', 'QC_after_trimming', label='FASTQC')

    graph.node('reference_genome', color=default_color)
    graph.edge('reference_genome', 'reference_genome', label="extract introns, \n STAR indexing")

    graph.node('BAM', color=default_color)
    graph.edge('FASTQ_trimmed', 'BAM', label='Alignment')
    graph.edge('reference_genome', 'BAM')

    graph.node('BAM_before_deduplication', color=umi_color)
    graph.edge('FASTQ_trimmed', 'BAM_before_deduplication', label='Alignment', color=umi_color)
    graph.edge('reference_genome', 'BAM_before_deduplication', color=umi_color)
    graph.edge('BAM_before_deduplication', 'BAM', label='Deduplication', color=umi_color)

    graph.node('feature_counts', color=default_color)
    graph.edge('BAM', 'feature_counts', label='featureCounts')
    graph.edge('reference_genome', 'feature_counts')

    graph.node('feature_counts_aggregated', color=default_color)
    graph.edge('feature_counts', 'feature_counts_aggregated', label='Aggregation')

    graph.node('DEX_analysis', color=default_color)
    graph.edge('feature_counts_aggregated', 'DEX_analysis')

    graph.node('coverage', color=default_color)
    graph.edge('BAM', 'coverage', label='Compute coverage')

    graph.node('intron_slopes', color=default_color)
    graph.edge('coverage', 'intron_slopes', label='Estimate slopes')
    graph.edge('reference_genome', 'intron_slopes')

    graph.node('selected_intron_slopes', color=default_color)
    graph.edge('intron_slopes', 'selected_intron_slopes', label='Select intron slopes')

    graph.edge('BAM', 'selected_intron_slopes', label='SJ info')
    graph.edge('BAM_before_deduplication', 'selected_intron_slopes', color=umi_color)

    graph.node('strandedness_info', color=default_color)
    graph.edge('BAM', 'strandedness_info', label='Infer strandedness')

    graph.view()

    (graph_for_render,) = pydot.graph_from_dot_file('dag.gv')
    graph_for_render.write_png('dag.png')
