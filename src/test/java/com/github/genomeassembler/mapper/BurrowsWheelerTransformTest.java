package com.github.genomeassembler.mapper;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class BurrowsWheelerTransformTest {

    @Test
    public void testTransform() {
        String test1 = "banana";
        String test2 = "panamabananas";
        String test3 = "mississippi";
        String test4 = "ACTACG";
        String expected1 = "annb$aa";
        String expected2 = "smnpbnnaaaaa$a";
        String expected3 = "ipssm$pissii";
        String expected4 = "GT$AACC";
        assertEquals(expected1, BurrowsWheelerTransform.transform(test1));
        assertEquals(expected2, BurrowsWheelerTransform.transform(test2));
        assertEquals(expected3, BurrowsWheelerTransform.transform(test3));
        assertEquals(expected4, BurrowsWheelerTransform.transform(test4));
    }

    @Test
    public void testInvert() {
        String test1 = "annb$aa";
        String test2 = "smnpbnnaaaaa$a";
        String test3 = "ipssm$pissii";
        String test4 = "GT$AACC";
        String expected1 = "banana";
        String expected2 = "panamabananas";
        String expected3 = "mississippi";
        String expected4 = "ACTACG";
        assertEquals(expected1, BurrowsWheelerTransform.invert(test1));
        assertEquals(expected2, BurrowsWheelerTransform.invert(test2));
        assertEquals(expected3, BurrowsWheelerTransform.invert(test3));
        assertEquals(expected4, BurrowsWheelerTransform.invert(test4));
    }
}
