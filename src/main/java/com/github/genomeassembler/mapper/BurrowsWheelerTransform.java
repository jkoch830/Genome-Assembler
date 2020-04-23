package com.github.genomeassembler.mapper;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;


/**
 * Contains static method for Burrows-Wheeler transformation and suffix
 * array construction
 */
public class BurrowsWheelerTransform {

    /**
     * Performs the Burrows Wheeler transform
     * @param genome The string being transformed
     *               - Requires that '$' is not present
     * @param suffixArray The buffer in which the suffix array will be placed
     *                    - Requires that size of suffix array is
     *                      is len(genome) + 1
     * @return The transformed string
     */
    public static String transform(String genome, int[] suffixArray) {
        // Initialize suffix array
        for (int i = 1; i < suffixArray.length; i++) {
            suffixArray[i - 1] = i;
        }
        suffixArray[suffixArray.length - 1] = 0;

        StringBuilder transformed = new StringBuilder(genome);
        transformed.append("$");
        int n = transformed.length(), p, r;
        char c;
        for (int s = n - 2; s >= 0; s--) {
            c = transformed.charAt(s);
            p = transformed.indexOf("$");
            r = s;
            for (int j = s + 1; j < n; j++) {
                char current = transformed.charAt(j);
                if (current < c || (j <= p && current == c)) {
                    r++;
                }
            }
            transformed.replace(p, p + 1, String.valueOf(c));
            transformed.deleteCharAt(s);
            transformed.insert(r, "$");
            swap(suffixArray, s, p);
            for (int k = s; k < r; k++) {
                swap(suffixArray, k, k + 1);
            }
        }
        return transformed.toString();
    }

    /**
     * Inverts a Burrows Wheeler transformed string back to its original state
     * @param transformed The transformed string
     * @return The original string
     */
    public static String invert(String transformed) {
        StringBuilder inverted = new StringBuilder("$");
        char[] tempArray = transformed.toCharArray();
        Arrays.sort(tempArray);
        String transformedSorted = String.copyValueOf(tempArray);
        tempArray = null;    // Garbage collect

        int lettersAdded = 0;
        int i = 0;
        while (lettersAdded < transformed.length() - 1) {
            char currChar = transformed.charAt(i);
            inverted.append(currChar);
            lettersAdded++;
            int numOccurrences = findNumOccurrences(transformed, currChar, i);
            i = findOccurrenceIndex(transformedSorted, currChar, numOccurrences);
        }
        inverted.deleteCharAt(inverted.indexOf("$"));
        inverted.reverse();
        return inverted.toString();
    }

    /**
     * Finds the number of occurrences of a character up to a certain index
     * @param string The string being searched
     * @param c The character whose occurrences are being recorded
     * @param end The range of searching (inclusive)
     * @return The number of occurrences of character c in string up to index end
     */
    private static int findNumOccurrences(String string, char c, int end) {
        int count = 0;
        for (int i = 0; i <= end; i++) {
            if (string.charAt(i) == c) {
                count++;
            }
        }
        return count;
    }

    /**
     * Gets index of string after a specific character has appeared a specific
     * amount of times
     * @param string The string
     * @param c The character being searched
     * @param occurrences The number of occurrences of c before returning
     * @return The index after 'c' has occurred 'occurrences' times
     */
    private static int findOccurrenceIndex(String string, char c, int occurrences) {
        int count = 0;
        for (int i = 0; i < string.length(); i++) {
            if (string.charAt(i) == c) {
                count++;
                if (count == occurrences) {
                    return i;
                }
            }
        }
        return -1;
    }

    private static void swap(int[] A, int i, int j) {
        int temp = A[i];
        A[i] = A[j];
        A[j] = temp;
    }


}
