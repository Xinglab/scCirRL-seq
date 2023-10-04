import sys
# import pyabpoa as pa
import edlib as ed
from collections import defaultdict as dd
import multiprocessing as mp

from . import seq_utils as su
from . import utils as ut

NULL_node_id = -1


class UMIDirectionalNetworkNode:
    def __init__(self, umi, reads, trans, node_id):
        self.umi = umi
        self.reads = reads
        self.node_id = node_id
        self.son_ids = []
        self.top_id = NULL_node_id
        self.cmpt_trans = trans
        # self.cDNA_seqs = []
        # self.poa_aln = pa.msa_aligner()
        # self.cons = ''

    def set_son_id(self, son_id):
        if son_id not in self.son_ids:
            self.son_ids.append(son_id)


# UMI deduplication criteria:
# 1. ED < T (2)
# 2. overlapping cmpt_trans
def search_for_umi_from_node(umi_graph_id_dict, node_id, umi_max_ed):
    node = umi_graph_id_dict[node_id]
    umi = node.umi
    min_ed = len(umi)
    min_intersection_trans = []
    min_node_id = -1
    # reads = node.reads
    for i in range(0, node_id):
        from_node = umi_graph_id_dict[i]
        # if len(from_node.reads) < 2 * len(node.reads)-1:
            # break
        from_umi = from_node.umi
        # 1. ed.align of two UMIs
        res = ed.align(umi, 'N'*umi_max_ed+from_umi+'N'*umi_max_ed, task="locations", mode="HW", k=umi_max_ed)
        if res['editDistance'] == -1:
            continue
        elif res['editDistance'] >= min_ed:
            continue

        # 2. overlapping compatible transcript set
        from_trans = from_node.cmpt_trans
        trans = node.cmpt_trans
        if from_trans == ['NA'] or trans == ['NA']:
            continue
        intersection_trans = list(set(from_trans).intersection(set(trans)))
        # print('{} {}'.format(from_trans, trans))
        if intersection_trans:
            min_ed = res['editDistance']
            min_intersection_trans = intersection_trans
            min_node_id = i
            # from_node.cmpt_trans = list(intersection_trans)
            # node.cmpt_trans = list(intersection_trans)
            # return i
    if min_node_id != -1:
        node.cmpt_trans = min_intersection_trans
        umi_graph_id_dict[min_node_id].cmpt_trans = list(min_intersection_trans)
        return min_node_id
    else:
        return NULL_node_id


# return: umi: {read_set: {trans_set}}
# allow same UMI but different trans
def get_umi_trans(umi, umi_reads, all_read_to_trans):
    read_to_trans = dict()
    for read in umi_reads:
        if read not in all_read_to_trans:
            read_to_trans[read] = ['NA']  # not listed in ESPRESSO, set trans as NA
        else:
            read_to_trans[read] = all_read_to_trans[read]
    if len(read_to_trans) == 0:
        return {}
    read_to_trans = dict(sorted(read_to_trans.items(), key=lambda d: (len(d[1]))))
    all_reads = list(read_to_trans.keys())
    read0 = all_reads[0]
    intersection_trans = set(read_to_trans[read0])
    for read in list(read_to_trans.keys())[1:]:
        trans = set(read_to_trans[read])
        intersection_trans = intersection_trans.intersection(trans)

    if not intersection_trans:  # empty trans set
        i = 1
        for read in read_to_trans:
            sys.stderr.write('Same BC/UMI & different isoforms: {}\t{} ({})\t{}\n'.format(read, umi, i, ','.join(read_to_trans[read])))
            i += 1
        return read_to_trans
    else:
        return {','.join(all_reads): intersection_trans}


def build_umi_network_graph(umi_to_reads, perfect_umi, read_to_trans, umi_max_ed):
    umi_graph_id_dict = dict()
    # sort by read count and sort by if umi/read is perferct_bc_umi
    umi_to_reads = dict(sorted(umi_to_reads.items(), key=lambda d: (-len(d[1]), d[0] not in perfect_umi)))
    node_id = 0
    for umi in umi_to_reads:
        umi_reads = umi_to_reads[umi]
        # same bc & UMI, but different transcript (ESPRESSO)
        umi_read_to_trans = get_umi_trans(umi, umi_reads, read_to_trans)
        if not umi_read_to_trans:
            continue
        for all_reads, trans in umi_read_to_trans.items():
            reads = all_reads.rsplit(',')
            node = UMIDirectionalNetworkNode(umi, reads, trans, node_id)

            umi_graph_id_dict[node_id] = node
            # search for from_node
            from_id = search_for_umi_from_node(umi_graph_id_dict, node_id, umi_max_ed)
            if from_id != NULL_node_id:
                from_node = umi_graph_id_dict[from_id]
                node.top_id = from_node.top_id  # XXX top_id: perfect_bc_umi with a higher read count XXX
                top_node = umi_graph_id_dict[node.top_id]
                top_node.cmpt_trans = from_node.cmpt_trans
                top_node.set_son_id(node_id)
            else:
                node.top_id = node_id
            node_id += 1
    return umi_graph_id_dict


def get_cluster_reads(top_node, umi_network_graph):
    tot_read = dict()
    for read in top_node.reads:
        tot_read[read] = 1
    for son_id in top_node.son_ids:
        son_node = umi_network_graph[son_id]
        for read in son_node.reads:
            tot_read[read] = 1
    return list(tot_read.keys())


# output NA for reads without compatible isoforms
def umi_clustering(read_to_trans, bu_res, umi_max_ed):
    if len(bu_res) == 0:
        return []
    umi_cluster_res = []
    bc_to_umi_reads = dd(lambda: dd(lambda: []))
    perfect_bc_to_umi = dd(lambda: [])
    # assume that any two reads with same BC & UMI always come from same transcript XXX
    # ESPRESSO/ISOQUANT may assign them with different transcripts
    for qname, bc, umi, is_perfect in bu_res:
        bc_to_umi_reads[bc][umi].append(qname)
        if is_perfect:
            perfect_bc_to_umi[bc].append(umi)
    for bc in bc_to_umi_reads:
        # if bc == 'TTGTGGATCATTCACT':
            # print('ok')
        umi_network_graph = build_umi_network_graph(bc_to_umi_reads[bc], perfect_bc_to_umi[bc], read_to_trans, umi_max_ed)
        for node_id, node in umi_network_graph.items():
            if node.top_id == NULL_node_id:
                ut.err_fatal_format_time(__name__, 'NULL top id: {}\n'.format(node.umi))
            top_id = node.top_id
            if top_id == node_id:  # only output top_id node
                umi = umi_network_graph[top_id].umi
                cmpt_trans = node.cmpt_trans
                reads = get_cluster_reads(umi_network_graph[top_id], umi_network_graph)
                umi_cluster_res.append((bc, umi, reads, cmpt_trans))
    return umi_cluster_res
