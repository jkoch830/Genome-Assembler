package com.github.genomeassembler;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
        assertEquals(expectedMappedReads, assembler.getMappedReads().keySet());
        assertEquals(Collections.singletonList(0), assembler.getMappedReads().get("ACTTG"));
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

        List<String> expectedUnmappedReads1 = new ArrayList<>(Collections.singletonList(
                "CGTAAGT"
        ));
        Set<String> expectedMappedReads1 = new HashSet<>(Arrays.asList(
                "ACGT", "GTACGTA", "AATT", "TTCC", "CCGG",
                "ACGTACGTAATTCCGG", "TACGTAA"
        ));
        assertEquals(8, assembler.mapReads(0));
        assertEquals(expectedUnmappedReads1, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads1, assembler.getMappedReads().keySet());

    }

    @Test
    public void testExactContigGenerationNoOverlap() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "CTAGATCGAT", "ATCGATCAGTCAC", "TATTA", "CTTAA"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        assembler.mapReads(0);
        assertEquals(2, assembler.formContigs());
        assertTrue(assembler.getMappedReads().isEmpty());
        Map<String, Integer> mappedContigs = assembler.getMappedContigSets().get(0);
        assertTrue(mappedContigs.containsKey("CTAGATCGATCAGTCACTATTA"));
    }

    @Test
    public void testExactContigGenerationOneOverlap() {
        String genome = "ACTAGATCGATCAGTCACTATTACCCTTAA";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "CTAGATCGAT", "ATCGATCAGTCAC", "TATTA", "CTTAA"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        AssemblerParameters param = new AssemblerParameters.Builder().requiredContigOverlap(1).build();
        assembler.setAssemblerParameters(param);
        assembler.mapReads(0);
        assertEquals(3, assembler.formContigs());
        Map<String, Integer> mappedContigs = assembler.getMappedContigSets().get(0);
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
        AssemblerParameters param = new AssemblerParameters.Builder().requiredContigOverlap(4).build();
        assembler.setAssemblerParameters(param);
        assembler.mapReads(0);
        assertEquals(4, assembler.formContigs());
        Map<String, Integer> mappedContigs = assembler.getMappedContigSets().get(0);
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
        AssemblerParameters param = new AssemblerParameters.Builder().requiredContigOverlap(1).build();
        assembler.setAssemblerParameters(param);
        assembler.mapReads(0);
        assertEquals(1, assembler.formContigs());
        Map<String, Integer> mappedContigs = assembler.getMappedContigSets().get(0);
        assertTrue(mappedContigs.containsKey(genome));
    }

    @Test
    public void testTwoErrorContigGenerationOverlap() {
        String genome = "TTGACTGAGCTGC";
        List<String> reads = new ArrayList<>(Arrays.asList(
                "ACCTA", "GTGAGCT"
        ));
        GenomeAssembler assembler = new GenomeAssembler(genome, reads);
        AssemblerParameters param = new AssemblerParameters.Builder().requiredContigOverlap(1).build();
        assembler.setAssemblerParameters(param);
        assembler.mapReads(2);
        assertEquals(1, assembler.formContigs());
        Map<String, Integer> mappedContigs = assembler.getMappedContigSets().get(0);
        int index = mappedContigs.get("ACCTAGCTA");
        assertEquals(3, index);
    }
}
