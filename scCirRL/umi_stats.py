"""Shared output schema and summaries for distinct raw UMIs per cluster."""

BC_UMI_FIELDS = [
    '#CellBarcode', 'UMI', 'ReadCount', 'TranscriptID', 'GeneID', 'GeneName',
    'ReadNames', 'TransCate', 'SpliceTag', 'RawUMICount',
]
BC_UMI_HEADER = '\t'.join(BC_UMI_FIELDS) + '\n'


def format_umi_clustering_stats(raw_umi_count_histogram, missing_clusters=0):
    """Format a histogram mapping distinct raw UMI counts to cluster counts."""
    counted_clusters = sum(raw_umi_count_histogram.values())
    raw_umis = sum(k * n for k, n in raw_umi_count_histogram.items())
    merged_clusters = sum(n for k, n in raw_umi_count_histogram.items() if k > 1)
    mean = raw_umis / counted_clusters if counted_clusters else 0.0
    merged_percent = 100.0 * merged_clusters / counted_clusters if counted_clusters else 0.0
    lines = [
        'UMI clustering statistics:',
        'Total UMI clusters                : {}'.format(counted_clusters + missing_clusters),
        'Clusters with raw UMI counts      : {}'.format(counted_clusters),
        'Clusters missing raw UMI counts   : {}'.format(missing_clusters),
        'Clusters with >=2 distinct UMIs   : {} ({:.2f}% of counted clusters)'.format(
            merged_clusters, merged_percent),
        'Raw UMIs summed within clusters   : {}'.format(raw_umis),
        'Additional raw UMIs merged        : {}'.format(raw_umis - counted_clusters),
        'Mean raw UMIs per counted cluster : {:.4f}'.format(mean),
        'Max raw UMIs per counted cluster  : {}'.format(max(raw_umi_count_histogram, default=0)),
        'Distinct raw UMIs per cluster     : cluster count (% of counted clusters)',
    ]
    for raw_count, cluster_count in sorted(raw_umi_count_histogram.items()):
        percent = 100.0 * cluster_count / counted_clusters if counted_clusters else 0.0
        lines.append('{:>32} : {} ({:.2f}%)'.format(raw_count, cluster_count, percent))
    return '\n'.join(lines)
