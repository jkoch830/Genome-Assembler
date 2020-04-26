package com.github.genomeassembler.mapper;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import static org.junit.Assert.assertArrayEquals;
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
        int[] suffixArray1 = new int[test1.length() + 1];
        int[] suffixArray2 = new int[test2.length() + 1];
        int[] suffixArray3 = new int[test3.length() + 1];
        int[] suffixArray4 = new int[test4.length() + 1];
        int[] expectedArray1 = getSuffixArray(test1);
        int[] expectedArray2 = getSuffixArray(test2);
        int[] expectedArray3 = getSuffixArray(test3);
        int[] expectedArray4 = getSuffixArray(test4);

        assertEquals(expected1, BurrowsWheelerTransform.transform(test1, suffixArray1));
        assertEquals(expected2, BurrowsWheelerTransform.transform(test2, suffixArray2));
        assertEquals(expected3, BurrowsWheelerTransform.transform(test3, suffixArray3));
        assertEquals(expected4, BurrowsWheelerTransform.transform(test4, suffixArray4));

        assertArrayEquals(expectedArray1, suffixArray1);
        assertArrayEquals(expectedArray2, suffixArray2);
        assertArrayEquals(expectedArray3, suffixArray3);
        assertArrayEquals(expectedArray4, suffixArray4);
    }

    @Test
    public void testInvert() {
        String test1 = "annb$aa";
        String test2 = "smnpbnnaaaaa$a";
        String test3 = "ipssm$pissii";
        String test4 = "GT$AACC";
        String test5 = "AAAAAAACCCCCCTTTTTTGGGGGGCATCGTAGCTAGCTGCTA";
        String expected1 = "banana";
        String expected2 = "panamabananas";
        String expected3 = "mississippi";
        String expected4 = "ACTACG";
        int[] suffixArray1 = new int[test5.length() + 1];
        assertEquals(expected1, BurrowsWheelerTransform.invert(test1));
        assertEquals(expected2, BurrowsWheelerTransform.invert(test2));
        assertEquals(expected3, BurrowsWheelerTransform.invert(test3));
        assertEquals(expected4, BurrowsWheelerTransform.invert(test4));

        assertEquals(test5, BurrowsWheelerTransform.invert(
                BurrowsWheelerTransform.transform(test5, suffixArray1)));
        assertArrayEquals(getSuffixArray(test5), suffixArray1);
    }

    /**
     * Naive (inefficient) correct implementation of suffix array construction
     * @param s The string
     * @return The suffix array
     */
    private static int[] getSuffixArray(String s) {
        s += "$";
        List<String> suffixes = new ArrayList<>();
        Map<String, Integer> mapping = new HashMap<>();
        for (int i = 0; i < s.length(); i++) {
            suffixes.add(s.substring(i));
            mapping.put(s.substring(i), i);
        }
        Collections.sort(suffixes);
        int[] suffixArray = new int[s.length()];
        for (int i = 0; i < suffixes.size(); i++) {
            suffixArray[i] = mapping.get(suffixes.get(i));
        }
        return suffixArray;
    }


}
