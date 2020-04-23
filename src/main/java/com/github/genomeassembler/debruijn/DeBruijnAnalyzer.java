package com.github.genomeassembler.debruijn;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Contains methods that perform graph analysis algorithms
 */
public class DeBruijnAnalyzer {

    private static boolean isOneInOneOut(DeBruijnGraph graph, String node) {
        return graph.getInDegree(node) == graph.getOutDegree(node) &&
                graph.getInDegree(node) == 1;
    }


    /**
     * Generates all contigs in the graph by finding all maximal
     * non-branching paths in the graph
     * @param graph The graph being searched
     * @return A list of all the contigs
     */
    public static List<String> contigGeneration(DeBruijnGraph graph) {
        List<String> contigs = new ArrayList<>();
        for (String node : graph) {
            if (!isOneInOneOut(graph, node)) {
                if (graph.getOutDegree(node) > 0) {
                    for (String neighbor : graph.getOutNeighbors(node)) {
                        StringBuilder contig = new StringBuilder(node);
                        while (isOneInOneOut(graph, neighbor)) {
                            contig.append(neighbor.charAt(neighbor.length() - 1));
                            neighbor = graph.getOutNeighbors(neighbor).get(0);
                        }
                        contig.append(neighbor.charAt(neighbor.length() - 1));
                        contigs.add(contig.toString());
                    }
                }
            }
        }
        List<List<String>> cycles = getCycles(graph);
        System.out.println(cycles);
        for (List<String> cycle : cycles) {
            String start = cycle.get(0);
            String curr = cycle.get(1);
            StringBuilder contig = new StringBuilder(start);
            while (!curr.equals(start)) {
                contig.append(curr.charAt(curr.length() - 1));
                curr = graph.getOutNeighbors(curr).get(0);
            }
            contigs.add(contig.toString());
        }
        return contigs;
    }



    /**
     * Retrieves all the isolated cycles in the graph
     * @param graph The graph being searched
     * @return A list of lists, each containing an isolated cycle
     */
    public static List<List<String>> getCycles(DeBruijnGraph graph) {
        List<List<String>> cycles = new ArrayList<>();
        Set<String> checkedNodes = new HashSet<>();
        for (String node : graph) {
            if (checkedNodes.contains(node)) { // Node was already checked
                continue;
            }
            checkedNodes.add(node);
            if (isOneInOneOut(graph, node)) { // Node is possibly in cycle
                List<String> newCycle = new ArrayList<>();
                newCycle.add(node);
                String currNode = graph.getOutNeighbors(node).get(0);
                // Traverse until it branches or until it reaches start
                while (isOneInOneOut(graph, currNode)) {
                    if (currNode.equals(node)) { // Cycle was made
                        cycles.add(newCycle);
                        break;
                    } else if (checkedNodes.contains(currNode)) { // No cycle
                        break;
                    }
                    newCycle.add(currNode);
                    checkedNodes.add(currNode);
                    currNode = graph.getOutNeighbors(currNode).get(0);
                }
            }
        }
        return cycles;
    }
}
