package com.github.genomeassembler;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
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

    @Test
    public void temp() {
        List<Map<Integer, String>> mapList = new ArrayList<>();
        Map<Integer, String> contigMap1 = new HashMap<>();
        Map<Integer, String> contigMap2 = new HashMap<>();
        Map<Integer, String> contigMap3 = new HashMap<>();
        Map<Integer, String> contigMap4 = new HashMap<>();
        contigMap1.put(0, "TTACTG"); contigMap1.put(9, "TGC");
        contigMap2.put(1, "TGCACT"); contigMap2.put(9, "TGC");
        contigMap3.put(2, "CCTGAG"); contigMap3.put(9, "CAT");
        contigMap4.put(4, "TGTG"); contigMap4.put(9, "CAT");
        mapList.add(contigMap1);
        mapList.add(contigMap2);
        mapList.add(contigMap3);
        mapList.add(contigMap4);

        assertEquals('T', getConsensusBase(0, mapList));
        assertEquals('T', getConsensusBase(1, mapList));
        assertEquals('A', getConsensusBase(2, mapList));
        assertEquals('C', getConsensusBase(3, mapList));
        assertEquals('T', getConsensusBase(4, mapList));
        assertEquals('G', getConsensusBase(5, mapList));
        assertEquals('T', getConsensusBase(6, mapList));
        assertEquals('G', getConsensusBase(7, mapList));
        assertEquals('0', getConsensusBase(8, mapList));
        assertEquals('T', getConsensusBase(9, mapList));
        assertEquals('G', getConsensusBase(10, mapList));
        assertEquals('C', getConsensusBase(11, mapList));


        Map<String, Integer> superContigs = new HashMap<>();
        int refGenomeLength = 12;
        StringBuilder superContig = new StringBuilder();
        for (int i = 0; i < refGenomeLength; i++) {
            char consensusBase = getConsensusBase(i, mapList);
            if (consensusBase == '0') { // No contig covers i
                if (superContig.length() != 0) { // don't save empty string
                    int startingIndex =  i - superContig.length();
                    superContigs.put(superContig.toString(), startingIndex);
                    superContig = new StringBuilder();
                }
            } else {
                superContig.append(consensusBase);
            }
        }

        if (superContig.length() != 0) {
            superContigs.put(superContig.toString(), refGenomeLength - superContig.length());
        }
        System.out.println(superContigs);
    }


    private static char getConsensusBase(int i, List<Map<Integer, String>> mapList) {
        List<Character> bases = new ArrayList<>();
        for (Map<Integer, String> map : mapList) {
            for (Map.Entry<Integer, String> entry : map.entrySet()) {
                int index = entry.getKey();
                String contig = entry.getValue();
                // Checks if contig covers i
                if (index <= i && i < index + contig.length()) {
                    char base = contig.charAt(i - index);
                    bases.add(base);
                    break;
                }
            }
        }
        if (bases.isEmpty()) { return '0'; } // if no contig covered i
        Map<Character, Integer> baseCount = new HashMap<>();
        for (char base : bases) { // Records frequencies
            if (!baseCount.containsKey(base)) {
                baseCount.put(base, 0);
            }
            baseCount.put(base, baseCount.get(base) + 1);
        }
        int maxCount = Collections.max(baseCount.values());
        for (char base : bases) { // breaks ties
            if (baseCount.get(base) == maxCount) {
                return base;
            }
        }
        return '0'; // Should never reach here
    }
}
