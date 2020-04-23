package com.github.genomeassembler.debuijn;

import com.github.genomeassembler.debruijn.BasicDeBruijnGraph;
import com.github.genomeassembler.debruijn.DeBruijnAnalyzer;
import com.github.genomeassembler.debruijn.DeBruijnGraph;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;

public class DeBruijnAnalyzerTest {

    @Test
    public void testFindCycles() {
        DeBruijnGraph graph = new BasicDeBruijnGraph();
        graph.addKmer("ACAT");
        graph.addKmer("CATC");
        graph.addKmer("ATCA");
        graph.addKmer("TCAC");
        graph.addKmer("CACA");

        graph.addKmer("AGAT");
        graph.addKmer("GATG");
        graph.addKmer("ATGA");
        graph.addKmer("TGAG");
        graph.addKmer("GAGA");

        List<List<String>> cycles = DeBruijnAnalyzer.getCycles(graph);
        assertEquals(2, cycles.size());
        assertEquals(5, cycles.get(0).size());
        assertEquals(5, cycles.get(1).size());

        graph.addKmer("AATG");
        cycles = DeBruijnAnalyzer.getCycles(graph);
        assertEquals(1, cycles.size());
        graph.addKmer("GATA");
        graph.addKmer("ATAA");
        graph.addKmer("TAAT");
        cycles = DeBruijnAnalyzer.getCycles(graph);
        assertEquals(1, cycles.size());


    }

    @Test
    public void testContigGeneration() {
        DeBruijnGraph graph = new BasicDeBruijnGraph();
        graph.addKmer("TAA");
        graph.addKmer("AAT");
        graph.addKmer("ATG");
        graph.addKmer("TGC");
        graph.addKmer("GCC");
        graph.addKmer("CCA");
        graph.addKmer("CAT");
        graph.addKmer("ATG");
        graph.addKmer("TGG");
        graph.addKmer("GGG");
        graph.addKmer("GGA");
        graph.addKmer("GAT");
        graph.addKmer("ATG");
        graph.addKmer("TGT");
        graph.addKmer("GTT");
        List<String> contigs = DeBruijnAnalyzer.contigGeneration(graph);
        assertEquals(9, contigs.size());
        List<String> expectedContigs = new ArrayList<>(Arrays.asList(
                "TAAT",
                "TGCCAT",
                "ATG",
                "ATG",
                "ATG",
                "GGAT",
                "TGG",
                "GGG",
                "TGTT"
        ));
        Collections.sort(contigs);
        Collections.sort(expectedContigs);
        assertEquals(expectedContigs, contigs);
        StringBuilder t = new StringBuilder("bananan");
        t.deleteCharAt(4);
        t.insert(6, "$");
        System.out.println(t);
    }

    @Test
    public void testContigGenerationLoops() {
        DeBruijnGraph graph = new BasicDeBruijnGraph();
        graph.addKmer("ACAT");
        graph.addKmer("CATC");
        graph.addKmer("ATCA");
        graph.addKmer("TCAC");
        graph.addKmer("CACA");

        graph.addKmer("AGAT");
        graph.addKmer("GATG");
        graph.addKmer("ATGA");
        graph.addKmer("TGAG");

        List<String> contigs = DeBruijnAnalyzer.contigGeneration(graph);
        assertEquals(2, contigs.size());
        assertTrue(contigs.contains("AGATGAG"));
        assertTrue(contigs.contains("ATCACAT"));
    }
}
