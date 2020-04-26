package com.github.genomeassembler.debruijn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Naive implementation of a de Bruijn graph. Nodes are stored as strings
 * and the graph is represented by an adjacency list. This is not suitable for
 * actual genome assembly due to the large memory requirements of strings.
 */
public class BasicDeBruijnGraph implements DeBruijnGraph {

    private Map<String, List<String>> graph;
    private Map<String, List<String>> inNeighbors;
    private int n = 0;
    private int m = 0;

    public BasicDeBruijnGraph() {
        graph = new HashMap<>();
        inNeighbors = new HashMap<>();
    }


    @Override
    public void addKmer(String kmer) {
        int k = kmer.length();
        String prefix = kmer.substring(0, k - 1);
        String suffix = kmer.substring(1);
        if (!graph.containsKey(prefix)) {
            n++;
            graph.put(prefix, new ArrayList<>());
            inNeighbors.put(prefix, new ArrayList<>());
        }
        graph.get(prefix).add(suffix); // Adds an edge from prefix to suffix
        if (!graph.containsKey(suffix)) { // Initializes suffix as node
            n++;
            graph.put(suffix, new ArrayList<>());
            inNeighbors.put(suffix, new ArrayList<>());
        }
        inNeighbors.get(suffix).add(prefix);
        m++;
    }

    @Override
    public List<String> getOutNeighbors(String node) {
        return (graph.containsKey(node)) ? graph.get(node) : new ArrayList<>();
    }

    @Override
    public List<String> getInNeighbors(String node) {
        return (inNeighbors.containsKey(node)) ? inNeighbors.get(node) : new ArrayList<>();
    }

    @Override
    public int getInDegree(String node) {
        return (inNeighbors.containsKey(node)) ? inNeighbors.get(node).size() : -1;
    }

    @Override
    public int getOutDegree(String node) {
        return (graph.containsKey(node)) ? graph.get(node).size() : -1;
    }

    @Override
    public int getNumNodes() {
        return n;
    }

    @Override
    public int getNumEdges() {
        return m;
    }

    @Override
    public Iterator<String> iterator() {
        return graph.keySet().iterator();
    }
}
