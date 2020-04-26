package com.github.genomeassembler.mapper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;


/**
 * The mapping utilizes the Burrows-Wheeler transformation. It stores the full
 * genome in an int array.
 * Each i-th int stores information of 3 items:
 *      - The first two bits store i-th char of the sorted BWT string
 *      - The 3rd and 4th bit store the i-th char of the BWT string
 *      - The other 28 bits store how many times the i-th char of the BWT string
 *              occurred within the BWT string up to index i (inclusive)
 * This does not support genomes longer than ~200,000,000 bp long
 * This costs ~4 * len(genome) bytes plus expensive Java overhead
 */
public class BWReadMapper implements ReadMapper {
    private final static int SORTED_BITS_OFFSET = 30;
    private final static int BWT_BITS_OFFSET = 28;
    private final static int CHAR_MASK = 0x3;
    private final static int MAX_OCCURRENCES = 0xFFFFFFF;
    private final static String BWT_FINISH_MSG = "Completed Burrows-Wheeler " +
            "transformation and suffix array construction";

    private int endMarkerLocation = -1;  // Marks where $ is
    private int numA;
    private int numC;
    private int numG;
    private int numT;
    private final int[] genome;
    private int[] suffixArray; // Stores every suffix


    /**
     * Processes the genome
     *      - Performs the Burrows Wheeler transformation
     *      - Counts number of A's, C's, T's, and G's within the genome
     *      - Stores information within an int[] as long as the genome + 1
     * @param genomeSequence The genome sequence
     */
    public BWReadMapper(String genomeSequence) {
        this.numA = 0;
        this.numC = 0;
        this.numT = 0;
        this.numG = 0;
        this.genome = new int[genomeSequence.length() + 1]; // plus a $
        //this.suffixArray = new int[genomeSequence.length() + 1]; // plus a $
        //String bwt = BurrowsWheelerTransform.transform(genomeSequence, this.suffixArray);
        String bwt = quickSetup();
        System.out.println(BWT_FINISH_MSG);
        char[] tempArray = bwt.toCharArray();
        Arrays.sort(tempArray);
        String sorted = String.copyValueOf(tempArray);
        tempArray = null;
        char bwtChar, sortedChar;
        int packed, numOccurrences;
        for (int i = 0; i < genome.length; i++) {
            bwtChar = bwt.charAt(i);
            sortedChar = sorted.charAt(i);
            switch(bwtChar) {
                case 'A':
                    this.numA++;
                    numOccurrences = this.numA;
                    break;
                case 'C':
                    this.numC++;
                    numOccurrences = this.numC;
                    break;
                case 'G':
                    this.numG++;
                    numOccurrences = this.numG;
                    break;
                case 'T':
                    this.numT++;
                    numOccurrences = this.numT;
                    break;
                case '$':
                    this.endMarkerLocation = i;
                    numOccurrences = 1;
                    break;
                default:
                    throw new IllegalArgumentException("Char not recognized : " + bwtChar);
            }
            packed = packInformation(sortedChar, bwtChar, numOccurrences);
            this.genome[i] = packed;
        }
    }

    @Override
    public List<Integer> mapRead(String read, int mismatches) {
        char start = read.charAt(read.length() - 1);  // search backwards
        int rangeStart, rangeEnd; // Indices to check (end-exclusive)
        if (mismatches == 0) {
            if (start == 'A') {
                rangeStart = 1;
                rangeEnd = rangeStart + this.numA;
            } else if (start == 'C') {
                rangeStart = this.numA + 1;
                rangeEnd = rangeStart + this.numC;
            } else if (start == 'G') {
                rangeStart = this.numA + this.numC + 1;
                rangeEnd = rangeStart + this.numG;
            } else if (start == 'T') {
                rangeStart = this.numA + this.numC + this.numG + 1;
                rangeEnd = rangeStart + this.numT;
            } else {
                throw new IllegalArgumentException("Illegal read: " + read);
            }
        } else { // Can afford at least one mismatch
            rangeStart = 1; // Don't start at '$'
            rangeEnd = genome.length;
        }

        List<Integer> validMappings = new ArrayList<>();
        for (int i = rangeStart; i < rangeEnd; i++) {
            int endIndex;
            if (getSortedChar(genome[i]) != start) { // already one mismatch
                endIndex = tryMapping(read, mismatches - 1, i);
            } else {
                endIndex = tryMapping(read, mismatches, i);
            }

            if (endIndex != -1) {
                validMappings.add(this.suffixArray[endIndex]);
            }
        }
        return validMappings;
    }


    private int tryMapping(String read, int mismatches, int start) {
        char[] bases = read.toCharArray();
        char base, bwtChar;
        int index = start;
        int numMismatches = 0, numOccurrences, packed;
        for (int i = bases.length - 2; i >= 0; i--) { // Start at 2nd to last
            if (index == 0) { // Reached end
                return -1;
            }
            packed = this.genome[index];
            numOccurrences = getNumOccurrences(packed);
            if (index == this.endMarkerLocation) {
                bwtChar = '$';
            } else {
                bwtChar = getBwtChar(packed);
            }
            base = bases[i];
            if (base != bwtChar) {
                numMismatches++;
                if (numMismatches > mismatches) { // Too many mismatches
                    return -1;
                }
            }
            index = getSortedOccurrenceIndex(bwtChar, numOccurrences);
        }
        if (index == 0) { // Read can't map to end of string
            return -1;
        } else {
            return index;
        }
    }

    /**
     * Tries mapping multiple reads to a genome
     * This method removes any mapped read from the original list of reads
     * @param reads The list of reads being mapped
     * @param mismatches The number of tolerant mismatches
     * @return Mappings of a reads to a list of starting positions in genome
     */
    public Map<String, List<Integer>> mapReads(List<String> reads, int mismatches) {
        Map<String, List<Integer>> mappings = new HashMap<>();
        for (String read : reads) {
            if (!mappings.containsKey(read)) {
                mappings.put(read, mapRead(read, mismatches));
            }
        }
        return mappings;
    }

    @Override
    public String getGenomeSequence() {
        StringBuilder originalGenome = new StringBuilder();
        for (int i = 0; i < genome.length; i++) {
            if (i == this.endMarkerLocation) {
                originalGenome.append('$');
            } else {
                originalGenome.append(getBwtChar(genome[i]));
            }
        }
        return BurrowsWheelerTransform.invert(originalGenome.toString());
    }

    /**
     * Returns a list of the counts of each base in the genome
     * @return The list of counts
     */
    public int[] getBaseCounts() {
        return new int[] {numA, numC, numG, numT};
    }

    /**
     * Packs information into an integer
     * @param sortedChar The char in the sorted BWT string
     * @param bwtChar The char in the BWT string
     * @param numOccurrences The number of times the char occurred in the BWT string
     * @return An int that's packed with this information
     */
    public static int packInformation(char sortedChar, char bwtChar, int numOccurrences) {
        if (numOccurrences > MAX_OCCURRENCES) {
            System.out.println("GENOME TOO LARGE");
        }
        Map<Character, Integer> encoding = new HashMap<>();
        encoding.put('A', 0);
        encoding.put('$', 0);
        encoding.put('C', 1);
        encoding.put('G', 2);
        encoding.put('T', 3);
        if (!encoding.containsKey(sortedChar) || !encoding.containsKey(bwtChar)) {
            throw new IllegalArgumentException("Char not recognized");
        }
        int packed = 0;
        packed |= (encoding.get(sortedChar) << SORTED_BITS_OFFSET);
        packed |= (encoding.get(bwtChar) << BWT_BITS_OFFSET);
        packed |= numOccurrences;
        return packed;
    }

    /**
     * Returns the 'sorted' char within the packed int
     * @param packed The packed int
     * @return The 'sorted' char
     */
    public static char getSortedChar(int packed) {
        return maskBits(packed, SORTED_BITS_OFFSET);
    }

    /**
     * Returns the 'BWT' char within the packed int
     * @param packed The packed int
     * @return The 'BWT' char
     */
    public static char getBwtChar(int packed) {
        return maskBits(packed, BWT_BITS_OFFSET);
    }

    /**
     * Retursn the number of occurrences within the packed int
     * @param packed The packed int
     * @return The number of occurrences stored inside the packed int
     */
    public static int getNumOccurrences(int packed) {
        return packed & MAX_OCCURRENCES;
    }

    private static char maskBits(int packed, int offset) {
        int bits = (packed >> offset) & CHAR_MASK;
        if (bits == 0) {
            return 'A';
        } else if (bits == 1) {
            return 'C';
        } else if (bits == 2) {
            return 'G';
        } else {
            return 'T';
        }
    }

    private int getSortedOccurrenceIndex(char c, int occurrence) {
        assert(getSortedChar(genome[0]) == 'A');
        if (c == '$') { return 0; }
        int offset;
        if (c == 'A') { offset = 0; }
        else if (c == 'C') { offset = this.numA; }
        else if (c == 'G') { offset = this.numA + this.numC; }
        else if (c == 'T') { offset = this.numA + this.numC + this.numG; }
        else { throw new IllegalArgumentException("Invalid character: " + c); }
        return offset + occurrence;
    }

    public String quickSetup() {
        try {
            String path = "src/main/resources/M30BWT.txt";
            BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8);
            String bwt = br.readLine().replaceAll("\\s", "");
            String array = br.readLine().replaceAll(
                    "\\s", "").replaceAll(
                            "\\[", "").replaceAll(
                                    "]", "");
            String[] splitArray = array.split(",");
            this.suffixArray = Stream.of(splitArray).mapToInt(Integer::parseInt).toArray();
            return bwt;
        } catch (Exception e) {
            System.out.println("Error reading file: " + e);
            return "";
        }
    }



}
