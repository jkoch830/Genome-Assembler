package com.github.genomeassembler;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
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
        int count = assembler.mapExactReads();
        assertEquals(expectedUnmappedReads, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads, assembler.getExactMappedReads().keySet());
        assertEquals(Collections.singletonList(0), assembler.getExactMappedReads().get("ACTTG"));
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
        assertEquals(5, assembler.mapExactReads());
        assertEquals(expectedUnmappedReads1, assembler.getUnmappedReads());
        assertEquals(expectedMappedReads1, assembler.getExactMappedReads().keySet());


        List<String> expectedUnmappedReads2 = new ArrayList<>(Arrays.asList(
                 "CGTAAGT"
        ));
        Set<String> expectedMappedReads2 = new HashSet<>(Arrays.asList(
                "ACGT", "GTACGTA", "AATT", "TTCC", "CCGG",
                "CCGGAATTACGTACGT", "TTACGTA", "GGAA"
        ));
        assertEquals(3, assembler.mapExactReverseComplementaryReads());





    }


    @Test
    public void testContigGeneration() {
    }
}
