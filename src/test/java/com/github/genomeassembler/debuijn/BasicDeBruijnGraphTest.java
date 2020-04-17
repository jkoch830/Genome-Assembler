package com.github.genomeassembler.debuijn;

import com.github.genomeassembler.debruijn.BasicDeBruijnGraph;
import com.github.genomeassembler.debruijn.DeBruijnGraph;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class BasicDeBruijnGraphTest {
    private DeBruijnGraph graph;

    @Before
    public void setUp() {
        graph = new BasicDeBruijnGraph();
    }

    @Test
    public void testAddNode() {
        assertEquals(-1, graph.getInDegree("ACGT"));
        assertEquals(-1, graph.getOutDegree("ACGT"));
    }
}
