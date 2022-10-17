#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as f:
        for i in f:
            yield next(f).rstrip()
            next(f).rstrip()
            next(f).rstrip()


def cut_kmer(read, kmer_size):
    for i in range(len(read)):
        if i+kmer_size <= len(read):
            yield read[i:i+kmer_size]
        else:
            break


def build_kmer_dict(fastq_file, kmer_size):
    dico = {}
    for i in read_fastq(fastq_file):
        for n in cut_kmer(i, kmer_size):
            if n not in dico:
                dico[n] = 1
            else:
                dico[n] += 1
    return dico


def build_graph(kmer_dict):
    G = nx.DiGraph()
    key = list(kmer_dict.keys())

    for i in range(len(key)):
        if i + 1 <= len(key):
            value1 = key[i][:-1]
            value2 = key[i][1:]
            w = kmer_dict[key[i]]
            G.add_edge(value1, value2, weight=w)
        else:
            break

    return G

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    nodes = list(graph.nodes())
    starting_nodes = []
    for i in range(len(nodes)):
        if len(list(graph.predecessors(nodes[i]))) == 0:
            starting_nodes.append(nodes[i])

    return starting_nodes


def get_sink_nodes(graph):
    nodes = list(graph.nodes())
    ending_nodes = []
    for i in range(len(nodes)):
        if len(list(graph.successors(nodes[i]))) == 0:
            ending_nodes.append(nodes[i])
    
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs_list = []
    for i in range(len(starting_nodes)):
        for n in range(len(ending_nodes)):
            if nx.has_path(graph, starting_nodes[i], ending_nodes[n]):
                for path in nx.all_simple_paths(graph, starting_nodes[i], ending_nodes[n]):
                    tmp = path[0]
                    for r in range(1, len(path)):
                        tmp = tmp + path[r][-1]
                    contigs_list.append((tmp, len(tmp)))
    return contigs_list

def save_contigs(contigs_list, output_file):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)if __name__ == '__main__':
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)

def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


if __name__ == '__main__':
    file = '/home/sdv/m2bi/aly/Metagenomic/debruijn-tp/data/eva71_two_reads.fq'
    kmer_size = 6
    kmer_dico = {}
    kmer_dico = build_kmer_dict(file, kmer_size)
    graph = build_graph(kmer_dico)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    for path in nx.all_simple_paths(graph, starting_nodes[1], ending_nodes[0]):
        tmp = path[0]
        for i in range(len(path)):
            tmp = tmp + path[i][-1]



