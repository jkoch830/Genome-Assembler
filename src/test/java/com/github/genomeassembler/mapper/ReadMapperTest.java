package com.github.genomeassembler.mapper;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class ReadMapperTest {
    @Before
    public void setUp() {

    }

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

}
