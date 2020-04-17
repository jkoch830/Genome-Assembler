package com.github.genomeassembler.debuijn;

import com.github.genomeassembler.debruijn.BasicDeBruijnGraph;
import com.github.genomeassembler.debruijn.DeBruijnGraph;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static org.junit.Assert.*;

public class BasicDeBruijnGraphTest {
    private BasicDeBruijnGraph graph;

    @Before
    public void setUp() {
        graph = new BasicDeBruijnGraph();
    }

    @Test
    public void testInitialState() {
        assertEquals(-1, graph.getInDegree("ACGT"));
        assertEquals(-1, graph.getOutDegree("ACGT"));
        assertTrue(graph.getOutNeighbors("ACGT").isEmpty());
        assertTrue(graph.getInNeighbors("ACGT").isEmpty());
        Iterator<String> iterator = graph.iterator();
        assertFalse(iterator.hasNext());
    }

    @Test
    public void testAdd() {
        graph.addKmer("ACGT");
        assertEquals(-1, graph.getInDegree("ACGT"));
        assertEquals(-1, graph.getOutDegree("ACGT"));
        assertEquals(1, graph.getInDegree("CGT"));
        assertEquals(0, graph.getInDegree("ACG"));
        assertEquals(1, graph.getOutDegree("ACG"));
        assertEquals(0, graph.getOutDegree("CGT"));
        assertTrue(graph.getInNeighbors("ACGT").isEmpty());
        assertTrue(graph.getOutNeighbors("ACGT").isEmpty());
        assertTrue(graph.getInNeighbors("ACG").isEmpty());
        assertEquals("ACG", graph.getInNeighbors("CGT").get(0));
        assertTrue(graph.getOutNeighbors("CGT").isEmpty());
        assertEquals("CGT", graph.getOutNeighbors("ACG").get(0));
        Iterator<String> graphIterator = graph.iterator();
        int n = 0;
        Set<String> nodes = new HashSet<>(2);
        while (graphIterator.hasNext()) {
            n++;
            nodes.add(graphIterator.next());
        }
        assertEquals(2, n);
        Set<String> expectedNodes = new HashSet<>(Arrays.asList("ACG", "CGT"));
        assertEquals(expectedNodes, nodes);
    }
}
