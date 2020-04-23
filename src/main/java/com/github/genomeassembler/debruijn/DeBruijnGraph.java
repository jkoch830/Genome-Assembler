package com.github.genomeassembler.debruijn;

import java.lang.instrument.Instrumentation;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;



/**
 * We have ~8,589,000,000 bytes to work with
 * Each string is ((2 * len(string) + 38) and if not divisible by 8, round up to nearest 8) bytes
 * 150 length string is 344 bytes
 * 40 length string is 120 bytes
 * About 4,000,000 kmers, which is 1,376,000,000
 * The interface for a de Bruijn graph implementation where nodes are represent
 * a the prefixes and suffixes of kmers. This interface does
 * not use any Node objects as an attempt to reduce memory consumption.
 */
public interface DeBruijnGraph extends Iterable<String>{

    /**
     * Adds a kmer to the graph by creating a node for the prefix (if it's not
     * already in the graph) and creating a node for the suffix (if it's not
     * already in the graph) and by adding an edge between the two nodes.
     * The length of both of the prefix and the suffix is k-1
     */
    void addKmer(String kmer);


    /**
     * Finds all the neighbors that a node has a directed edge to
     * @param node The node
     * @return The list of neighbors, or an empty list if the node is absent
     */
    List<String> getOutNeighbors(String node);


    /**
     * Finds all the neighbors that have a directed edge to the given node
     * @param node The node
     * @return The list of neighbors, or an empty list if the node is absent
     */
    List<String> getInNeighbors(String node);


    /**
     * Retrieves the number of nodes that are directed toward a given kmer
     * @param node The kmer whose in-degree is being found
     * @return The in-degree of the kmer, or -1 if the node is not found
     */
    int getInDegree(String node);


    /**
     * Retrieves the number of nodes a given kmer is directed to
     * @param node The kmer whose out-degree is being found
     * @return The out-degree of the kmer, or -1 if the node is not found
     */
    int getOutDegree(String node);


    /**
     * Retrieves the number of nodes in the graph
     * @return The number of nodes
     */
    int getNumNodes();


    /**
     * Retrieves the number of edges in the graph
     * @return The number of edges
     */
    int getNumEdges();

}
