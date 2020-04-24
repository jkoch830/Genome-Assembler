package com.github.genomeassembler.mapper;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ReadMapperTest {

    @Test
    public void testPacking() {
        int packed1 = ReadMapper.packInformation('A', 'C', 5);
        int packed2 = ReadMapper.packInformation('T', 'G', 13232132);
        int packed3 = ReadMapper.packInformation('A', 'A', 0);
        assertEquals('A', ReadMapper.getSortedChar(packed1));
        assertEquals('T', ReadMapper.getSortedChar(packed2));
        assertEquals('A', ReadMapper.getSortedChar(packed3));
        assertEquals('C', ReadMapper.getBwtChar(packed1));
        assertEquals('G', ReadMapper.getBwtChar(packed2));
        assertEquals('A', ReadMapper.getBwtChar(packed3));
        assertEquals(5, ReadMapper.getNumOccurrences(packed1));
        assertEquals(13232132, ReadMapper.getNumOccurrences(packed2));
        assertEquals(0, ReadMapper.getNumOccurrences(packed3));
    }

    @Test
    public void testCounts() {
        ReadMapper mapper1 = new ReadMapper("ACATAGCTAGCTAGCA");
        ReadMapper mapper2 = new ReadMapper("AAACCC$$$");
        assertArrayEquals(new int[] {6, 4, 3, 3}, mapper1.getBaseCounts());
        assertArrayEquals(new int[] {3, 3, 0, 0}, mapper2.getBaseCounts());
    }

    @Test
    public void testGetGenomeSequence() {
        String genomeSequence1 = "ACGTTGCTAGCTAGCTAGCTAG";
        String genomeSequence2 = "AAAAAAACCCCCCTTTTTTGGGGGGCATCGTAGCTAGCTGCTA";
        String genomeSequence3 = "TTTTTTTTTTTTTTTT";
        String genomeSequence4 = "ACTGGTCAATCGGCACACGTAGCTAGCTATGGGGGGGGGGGG";
        ReadMapper mapper1 = new ReadMapper(genomeSequence1);
        ReadMapper mapper2 = new ReadMapper(genomeSequence2);
        ReadMapper mapper3 = new ReadMapper(genomeSequence3);
        ReadMapper mapper4 = new ReadMapper(genomeSequence4);
        assertEquals(genomeSequence1, mapper1.getGenomeSequence());
        assertEquals(genomeSequence2, mapper2.getGenomeSequence());
        assertEquals(genomeSequence3, mapper3.getGenomeSequence());
        assertEquals(genomeSequence4, mapper4.getGenomeSequence());


    }

    @Test
    public void testPerfectMapReads() {
        ReadMapper mapper1 = new ReadMapper("ACACTAGTCGATG");
        ReadMapper mapper2 = new ReadMapper("TTTTTCCCACATCATCATCA");
        List<Integer> mappings1 = mapper1.mapRead("AC", 0);
        List<Integer> mappings2 = mapper1.mapRead("A", 0);
        List<Integer> mappings3 = mapper1.mapRead("ACACT", 0);
        List<Integer> mappings4 = mapper1.mapRead("GA", 0);
        List<Integer> mappings5 = mapper1.mapRead("TAGT", 0);
        List<Integer> mappings6 = mapper1.mapRead("TT", 0);
        List<Integer> mappings7 = mapper2.mapRead("T", 0);
        List<Integer> mappings8 = mapper2.mapRead("TT", 0);
        List<Integer> mappings9 = mapper2.mapRead("TTCC", 0);
        List<Integer> mappings10 = mapper2.mapRead("CAT", 0);

        List<Integer> expected1 = new ArrayList<>(Arrays.asList(0, 2));
        List<Integer> expected2 = new ArrayList<>(Arrays.asList(0, 2, 5, 10));
        List<Integer> expected3 = new ArrayList<>(Collections.singletonList(0));
        List<Integer> expected4 = new ArrayList<>(Collections.singletonList(9));
        List<Integer> expected5 = new ArrayList<>(Collections.singletonList(4));
        List<Integer> expected7 = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 11, 14, 17));
        List<Integer> expected8 = new ArrayList<>(Arrays.asList(0, 1, 2, 3));
        List<Integer> expected9 = new ArrayList<>(Collections.singletonList(3));
        List<Integer> expected10 = new ArrayList<>(Arrays.asList(9, 12, 15));

        assertEquals(expected1, mappings1);
        assertEquals(expected2, mappings2);
        assertEquals(expected3, mappings3);
        assertEquals(expected4, mappings4);
        assertEquals(expected5, mappings5);
        assertTrue(mappings6.isEmpty());

        Collections.sort(mappings7);
        Collections.sort(mappings8);
        Collections.sort(mappings9);
        Collections.sort(mappings10);
        assertEquals(expected7, mappings7);
        assertEquals(expected8, mappings8);
        assertEquals(expected9, mappings9);
        assertEquals(expected10, mappings10);

    }

    @Test
    public void testImperfectMapReads() {
        ReadMapper mapper1 = new ReadMapper("ACTGCTTGT");
        ReadMapper mapper2 = new ReadMapper("AAAAATTTCTAGCTA");
        ReadMapper mapper3 = new ReadMapper("ACTGACTGGTCAG");
        ReadMapper mapper4 = new ReadMapper(
                "CATGCTGATCGTGATCGTAGCTAGTCGATCATGCTACTGGTCA"
        );

        List<Integer> mappings1 = mapper1.mapRead("AGT", 1);
        List<Integer> mappings2 = mapper1.mapRead("A", 2);
        List<Integer> mappings3 = mapper1.mapRead("ACTGCTAGT", 1);
        List<Integer> mappings4 = mapper1.mapRead("ACTACATAT", 3);
        List<Integer> mappings5 = mapper2.mapRead("AGT", 3);
        List<Integer> mappings6 = mapper2.mapRead("TTCCAA", 3);
        List<Integer> mappings7 = mapper3.mapRead("ACTGGCTAG", 2);
        List<Integer> mappings8 = mapper3.mapRead("TCA", 2);
        List<Integer> mappings9 = mapper3.mapRead("CTGAC", 4);
        List<Integer> mappings10 = mapper4.mapRead("CATGCTGATCG", 4);
        List<Integer> mappings11 = mapper4.mapRead("GATCGTAGCTAGTCGAT", 2);

        List<Integer> expected1 = new ArrayList<>(Arrays.asList(0, 6));
        List<Integer> expected2 = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8));
        List<Integer> expected3 = new ArrayList<>(Collections.singletonList(0));
        List<Integer> expected4 = new ArrayList<>(Collections.singletonList(0));
        List<Integer> expected5 = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12));
        List<Integer> expected6 = new ArrayList<>(Arrays.asList(5, 6, 9));
        List<Integer> expected7 = new ArrayList<>(Arrays.asList(0, 4));
        List<Integer> expected8 = new ArrayList<>(Arrays.asList(0, 2, 4, 6, 9));
        List<Integer> expected9 = new ArrayList<>(Arrays.asList(1, 5, 6, 8));
        List<Integer> expected10 = new ArrayList<>(Arrays.asList(0, 6, 29));
        List<Integer> expected11 = new ArrayList<>(Collections.singletonList(12));

        Collections.sort(mappings1);
        Collections.sort(mappings2);
        Collections.sort(mappings4);
        Collections.sort(mappings5);
        Collections.sort(mappings6);
        Collections.sort(mappings7);
        Collections.sort(mappings8);
        Collections.sort(mappings9);
        Collections.sort(mappings10);
        Collections.sort(mappings11);

        assertEquals(expected1, mappings1);
        assertEquals(expected2, mappings2);
        assertEquals(expected3, mappings3);
        assertEquals(expected4, mappings4);
        assertEquals(expected5, mappings5);
        assertEquals(expected6, mappings6);
        assertEquals(expected7, mappings7);
        assertEquals(expected8, mappings8);
        assertEquals(expected9, mappings9);
        assertEquals(expected10, mappings10);
        assertEquals(expected11, mappings11);


    }

    @Test
    public void testMapReads() {
        ReadMapper mapper1 = new ReadMapper("ACACTAGTCGATG");
        ReadMapper mapper2 = new ReadMapper("TTTTTTTTTTTT");
        List<Integer> mappings1 = mapper1.mapRead("CAC", 0);
        List<Integer> mappings2 = mapper1.mapRead("GATGA", 0);
        List<Integer> mappings3 = mapper1.mapRead("ACACTAGTCGATG", 0);
        List<Integer> mappings4 = mapper1.mapRead("ACACTAGTCGATGT", 0);
        List<Integer> mappings5 = mapper2.mapRead("TTTA", 1);

        List<Integer> expected1 = new ArrayList<>(Collections.singletonList(1));
        List<Integer> expected3 = new ArrayList<>(Collections.singletonList(0));
        List<Integer> expected5 = new ArrayList<>(Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8));

        Collections.sort(mappings5);

        assertEquals(expected1, mappings1);
        assertTrue(mappings2.isEmpty());
        assertEquals(expected3, mappings3);
        assertTrue(mappings4.isEmpty());
        assertEquals(expected5, mappings5);
    }

    @Test(expected=IllegalArgumentException.class)
    public void testIllegalCharacter() {
        ReadMapper mapper = new ReadMapper("ACATAGCTAGCRTAGCA");
    }


    private static void printArray(int[] array) {
        for (int n : array) {
            System.out.print(n + ", ");
        }
        System.out.println();
    }
}
