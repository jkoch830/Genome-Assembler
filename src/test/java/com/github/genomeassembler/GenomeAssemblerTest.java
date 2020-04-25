package com.github.genomeassembler;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
public class GenomeAssemblerTest {

    @Test
    public void testExactReadMapping() {
        String genome = "ACTTGCGTAGCTTGCTGATGT";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "ACTTG", "GCGTAGCTTGC", "GCG", "TGATGT",
                "CAGGATC", "CCATTGA", "GATG", "ATGT",
                "ACTTG"
        ));
        List<String> expectedUnmappedReads = new ArrayList<>(Arrays.asList(
                "CAGGATC", "CCATTGA"
        ));
        Set<String> expectedMappedReads = new HashSet<>(Arrays.asList(
                "ACTTG", "GCGTAGCTTGC", "GCG", "TGATGT", "GATG",
                "ATGT", "ACTTG"
                ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assertEquals(reads, assembler.getUnmappedReads());
        int count = assembler.mapReads(0);
        assertEquals(expectedUnmappedReads, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads, assembler.getMappedReads(0).keySet());
        assertEquals(Collections.singletonList(0), assembler.getMappedReads(0).get("ACTTG"));
        assertEquals(7, count);
    }

    @Test
    public void testReverseComplementaryReadMapping() {
        String genome = "ACGTACGTAATTCCGG";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "ACGT", "GTACGTA", "AATT", "TTCC", "CCGG",
                "GGAA", "TTACGTA", "CCGGAATTACGTACGT",
                "CGTAAGT"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);

        List<String> expectedUnmappedReads1 = new ArrayList<>(Arrays.asList(
                "GGAA", "TTACGTA", "CCGGAATTACGTACGT", "CGTAAGT"
        ));
        Set<String> expectedMappedReads1 = new HashSet<>(Arrays.asList(
                "ACGT", "GTACGTA", "AATT", "TTCC", "CCGG"
        ));
        assertEquals(5, assembler.mapReads(0));
        assertEquals(expectedUnmappedReads1, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads1, assembler.getMappedReads(0).keySet());

        // Map reverse complements
        List<String> expectedUnmappedReads2 = new ArrayList<>(Collections.singletonList(
                "CGTAAGT"
        ));
        Set<String> expectedMappedReads2 = new HashSet<>(Arrays.asList(
                "ACGT", "GTACGTA", "AATT", "TTCC", "CCGG",
                "ACGTACGTAATTCCGG", "TACGTAA"
        ));
        assertEquals(3, assembler.mapReverseComplementaryReads(0));
        assertEquals(expectedUnmappedReads2, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads2, assembler.getMappedReads(0).keySet());
    }

    @Test
    public void testExactContigGenerationNoOverlap() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "CTAGATCGAT", "ATCGATCAGTCAC", "TATTA", "CTTAA"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assembler.mapReads(0);
        assertEquals(2, assembler.formContigs(0, 0));
        assertTrue(assembler.getMappedReads(0).isEmpty());
        Map<String, Integer> mappedContigs = assembler.getMappedContigs(0);
        assertTrue(mappedContigs.containsKey("CTAGATCGATCAGTCACTATTA"));
    }

    @Test
    public void testExactContigGenerationOneOverlap() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "CTAGATCGAT", "ATCGATCAGTCAC", "TATTA", "CTTAA"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assembler.mapReads(0);
        assertEquals(3, assembler.formContigs(0, 1));
        Map<String, Integer> mappedContigs = assembler.getMappedContigs(0);
        assertTrue(mappedContigs.containsKey("CTAGATCGATCAGTCAC"));
        assertTrue(mappedContigs.containsKey("TATTA"));
        assertTrue(mappedContigs.containsKey("CTTAA"));
    }

    @Test
    public void testExactContigGenerationOverlap() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "CTAGATCGAT", "GATCAGTCAC", "TATTA", "CTTAA"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assembler.mapReads(0);
        assertEquals(4, assembler.formContigs(0, 4));
        Map<String, Integer> mappedContigs = assembler.getMappedContigs(0);
        assertTrue(mappedContigs.containsKey("CTAGATCGAT"));
        assertTrue(mappedContigs.containsKey("GATCAGTCAC"));
        assertTrue(mappedContigs.containsKey("TATTA"));
        assertTrue(mappedContigs.containsKey("CTTAA"));
    }

    @Test
    public void testExactContigGenerationFullCoverage() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "ACTAGATCGA", "CTA", "TCGATCA", "ATCAGTCA", "CACTATTAC",
                "ATTACCCTTA", "CCTTAA", "ATCAGTCAC"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assembler.mapReads(0);
        assertEquals(1, assembler.formContigs(0, 1));
        Map<String, Integer> mappedContigs = assembler.getMappedContigs(0);
        assertTrue(mappedContigs.containsKey(genome));
    }
}
