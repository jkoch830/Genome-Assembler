package com.github.genomeassembler.mapper;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


/**
 * Contains static method for Burrows-Wheeler transformation and suffix
 * array construction
 */
public class BurrowsWheelerTransform {

    /**
     * Performs the Burrows Wheeler transform
     * @param string The string being transformed
     *               - Requires that '$' is not present
     * @param suffixArray The buffer in which the suffix array will be placed
     *                    - Requires that size of suffix array is
     *                      is len(genome) + 1
     * @return The transformed string
     */
    public static String transform(String string, int[] suffixArray) {
        // Initialize suffix array
        for (int i = 1; i < suffixArray.length; i++) {
            suffixArray[i - 1] = i;
        }
        suffixArray[suffixArray.length - 1] = 0;
        string += "$";
        char[] transformed = string.toCharArray();
        int n = transformed.length, p = transformed.length - 1, r;
        char c, prevChar = transformed[n - 1];
        Map<Character, Integer> afterS = new HashMap<>();
        int percentUpdateIncrement = (n - 1 < 20) ? 1 : (n - 1) / 20;
        int percentComplete;
        for (int s = n - 2; s >= 0; s--) {
            if (n > 200 && (n - 2 - s) % percentUpdateIncrement == 0) {
                percentComplete = 5 * ((n - 2 - s) / percentUpdateIncrement);
                System.out.println(percentComplete + "% complete");
            }

            c = transformed[s];
            if (!afterS.containsKey(prevChar)) {
                afterS.put(prevChar, 0);
            }
            afterS.put(prevChar, afterS.get(prevChar) + 1);

            r = s;
            for (int j = s + 1; j < p; j++) {
                if (c == transformed[j]) {
                    r++;
                }
            }
            r += strictlySmaller(c, afterS);
            swapInts(suffixArray, s, p);
            swapChars(transformed, s, p);
            for (int k = s; k < r; k++) {
                swapInts(suffixArray, k, k + 1);
                swapChars(transformed, k, k + 1);
            }
            p = r;
            prevChar = c;
        }

        return String.valueOf(transformed);
    }


    private static int strictlySmaller(char c, Map<Character, Integer> map) {
        int count = 0;
        for (Map.Entry<Character, Integer> entry : map.entrySet()) {
            if (c > entry.getKey()) {
                count += entry.getValue();
            }
        }
        return count;
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

    private static void swapInts(int[] A, int i, int j) {
        int temp = A[i];
        A[i] = A[j];
        A[j] = temp;
    }

    private static void swapChars(char[] A, int i, int j) {
        char temp = A[i];
        A[i] = A[j];
        A[j] = temp;
    }


}
