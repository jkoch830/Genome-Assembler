package com.github.genomeassembler;

import com.github.genomeassembler.mapper.ReadMapper;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class GenomeAssembler {
    private final static String PROCESS_BEGIN_MSG = "Beginning processing of reference genome...";
    private final static String PROCESS_FINISH_MSG = "Finished processing of reference genome";
    private final static String MISMATCH_ERROR_MSG = "Number of mismatches must be between 0 and 5";
    private final static int REQUIRED_OVERLAP = 0;

    private final ReadMapper referenceGenomeReadMapper;

    private final List<String> unmappedReads;
    private final Map<String, List<Integer>> exactMappedReads;
    private final Map<String, List<Integer>> oneErrorMappedReads;
    private final Map<String, List<Integer>> twoErrorMappedReads;
    private final Map<String, List<Integer>> threeErrorMappedReads;
    private final Map<String, List<Integer>> fourErrorMappedReads;
    private final Map<String, List<Integer>> fiveErrorMappedReads;

    // Each correspond to a certain mismatch-tolerated mappedReads mapping
    private final Map<String, Integer> exactContigs;
    private final Map<String, Integer> oneErrorContigs;
    private final Map<String, Integer> twoErrorContigs;
    private final Map<String, Integer> threeErrorContigs;
    private final Map<String, Integer> fourErrorContigs;
    private final Map<String, Integer> fiveErrorContigs;


    public GenomeAssembler(String referenceGenomeSequence, List<String> reads) {
        // Process genome
        System.out.println("Length of genome: " + referenceGenomeSequence.length());
        System.out.println(PROCESS_BEGIN_MSG);
        this.referenceGenomeReadMapper = new ReadMapper(referenceGenomeSequence);
        System.out.println(PROCESS_FINISH_MSG);

        // Initialize read data structures
        this.unmappedReads = reads;
        this.exactMappedReads = new HashMap<>();
        this.oneErrorMappedReads = new HashMap<>();
        this.twoErrorMappedReads = new HashMap<>();
        this.threeErrorMappedReads = new HashMap<>();
        this.fourErrorMappedReads = new HashMap<>();
        this.fiveErrorMappedReads = new HashMap<>();

        // Initialize contig data structures
        this.exactContigs = new HashMap<>();
        this.oneErrorContigs = new HashMap<>();
        this.twoErrorContigs = new HashMap<>();
        this.threeErrorContigs = new HashMap<>();
        this.fourErrorContigs = new HashMap<>();
        this.fiveErrorContigs = new HashMap<>();

    }


    /**
     * Attempts to map every read with no tolerant mismatches.
     * If a read maps exactly, the read is removed from the unmappedReads list
     * to the mappedReads map
     * @param mismatches The number of tolerated mismatches when mapping
     * @return The number of reads (including duplicates) that were mapped
     *         to the reference genome
     */
    public int mapReads(int mismatches) {
        int count = 0;
        for (int i = this.unmappedReads.size() - 1; i >= 0; i--) {
            String read = this.unmappedReads.get(i);
            List<Integer> startingPositions =
                    this.referenceGenomeReadMapper.mapRead(read, mismatches);
            if (!startingPositions.isEmpty()) {
                this.unmappedReads.remove(i);
                exactMappedReads.put(read, startingPositions);
                count++;
            }
        }
        return count;
    }


    /**
     * Attempts to map remaining unmapped reads' reverse complements to the
     * reference genome with no tolerant mismatches.
     * If a read's reverse complement is mapped, the read is removed from
     * the unmapped reads, and its reverse complement is added to mapped
     * @param mismatches The number of tolerated mismatches when mapping
     * @return The number of reverse complementary reads (including duplicates)
     *         that were mapped to the reference genome
     */
    public int mapReverseComplementaryReads(int mismatches) {
        int count = 0;
        for (int i = this.unmappedReads.size() - 1; i >= 0; i--) {
            String reverseComplement = reverseComplement(this.unmappedReads.get(i));
            List<Integer> startingPositions =
                    this.referenceGenomeReadMapper.mapRead(reverseComplement, mismatches);
            if (!startingPositions.isEmpty()) {
                this.unmappedReads.remove(i);
                exactMappedReads.put(reverseComplement, startingPositions);
                count++;
            }
        }
        return count;
    }


    /**
     * Forms contigs from mapped reads that form a contiguous sequence.
     * The mapping containing the reads will be cleared because the resulting
     * contig map will contain reads that don't form a contig
     * @param mismatches The number of tolerated mismatches for the mapped reads
     *                   being converted to contigs.
     *                   Requires 0 <= mismatches <= 5
     * @param requiredOverlap The required amount of overlapping base pairs
     *                        two reads require in order to be combined into a
     *                        contig
     * @return The number of contigs formed
     */
    public int formContigs(int mismatches, int requiredOverlap) {
        Map<String, List<Integer>> mappedReads = getMappedReads(mismatches);
        if (mappedReads.isEmpty()) {
            return 0;
        }
        Map<Integer, String> mappedIndices = removeSmallStrings(flipMapping(mappedReads));
        mappedReads.clear(); // No need to store reads any longer

        Map<String, Integer> contigs = getMappedContigs(mismatches);

        Set<Integer> tempKeySet = mappedIndices.keySet();
        int start = Collections.min(tempKeySet);
        int end = Collections.max(tempKeySet);
        tempKeySet = null;

        String read, otherRead;
        int readLength;
        ContigBuilder newContig;
        for (int i = start; i <= end; i++) { // Loop through mapped indices
            if (mappedIndices.containsKey(i)) { // Found read at index
                read = mappedIndices.get(i);
                readLength = read.length();
                newContig = new ContigBuilder(read, i);
                // Tries finding other strings it overlaps with
                int j = i + 1;
                int validSearchBound = i + readLength - requiredOverlap, newBound;
                while (j <= validSearchBound) {
                    if (mappedIndices.containsKey(j)) { // overlaps with other read
                        otherRead = mappedIndices.get(j);
                        mappedIndices.remove(j); // remove other read from map
                        newContig.combineRead(otherRead, j);

                        newBound = j + otherRead.length() - requiredOverlap;
                        validSearchBound = Math.max(validSearchBound, newBound);
                    }
                    j++;
                }
                contigs.put(newContig.build(), i);
            }
        }
        return contigs.size();
    }


    /**
     * Retrieves the assembler's unmapped reads
     * @return A list containing all unmapped reads
     */
    public List<String> getUnmappedReads() {
        return this.unmappedReads;
    }


    /**
     * Retrieves the assembler's already mapped reads
     * @param mismatches The number of mismatches between the mapped reads and
     *                   the reference genome
     * @return A map containing each mapped read and their respective starting
     *      positions in the reference genome
     */
    public Map<String, List<Integer>> getMappedReads(int mismatches) {
        switch(mismatches) {
            case 0:
                return this.exactMappedReads;
            case 1:
                return this.oneErrorMappedReads;
            case 2:
                return this.twoErrorMappedReads;
            case 3:
                return this.threeErrorMappedReads;
            case 4:
                return this.fourErrorMappedReads;
            case 5:
                return this.fiveErrorMappedReads;
            default:
                throw new IllegalArgumentException(MISMATCH_ERROR_MSG);
        }
    }

    /**
     * Retrieve's the assembler's already formed contigs
     * @param mismatches The number of mismatches between the reads that formed
     *                   these contigs and the reference genome
     * @return A mapping between the contigs and their starting position in the
     *         genome
     */
    public Map<String, Integer> getMappedContigs(int mismatches) {
        switch(mismatches) {
            case 0:
                return this.exactContigs;
            case 1:
                return this.oneErrorContigs;
            case 2:
                return this.twoErrorContigs;
            case 3:
                return this.threeErrorContigs;
            case 4:
                return this.fourErrorContigs;
            case 5:
                return this.fiveErrorContigs;
            default:
                throw new IllegalArgumentException(MISMATCH_ERROR_MSG);
        }
    }

    private static String reverseComplement(String s) {
        StringBuilder complement = new StringBuilder();
        for (char c : s.toCharArray()) {
            if (c == 'A') {
                complement.insert(0, 'T');
            } else if (c == 'T') {
                complement.insert(0, 'A');
            } else if (c == 'C') {
                complement.insert(0, 'G');
            } else if (c == 'G') {
                complement.insert(0, 'C');
            } else {
                throw new IllegalArgumentException("Invalid character: " + c);
            }
        }
        return complement.toString();
    }

    private static Map<Integer, List<String>> flipMapping(Map<String,
            List<Integer>> mappedReads) {
        Map<Integer, List<String>> flippedMap = new HashMap<>();
        for (Map.Entry<String, List<Integer>> entry : mappedReads.entrySet()) {
            for (int startingPosition : entry.getValue()) {
                if (!flippedMap.containsKey(startingPosition)) {
                    flippedMap.put(startingPosition, new ArrayList<>());
                }
                flippedMap.get(startingPosition).add(entry.getKey());
            }
        }
        return flippedMap;
    }

    private static Map<Integer, String> removeSmallStrings(Map<Integer,
            List<String>> mappedIndices) {
        Map<Integer, String> simplified = new HashMap<>();
        for (Map.Entry<Integer, List<String>> entry : mappedIndices.entrySet()) {
            String longest = entry.getValue().get(0);
            int maxLength = longest.length();
            for (String s : entry.getValue()) { // Finds longest string
                if (s.length() > maxLength) {
                    maxLength = s.length();
                    longest = s;
                }
            }
            simplified.put(entry.getKey(), longest);
        }
        return simplified;
    }

    private static class ContigBuilder {
        private StringBuilder stringBuilder;
        private final int startIndex;

        private ContigBuilder(String read, int startIndex) {
            this.stringBuilder = new StringBuilder(read);
            this.startIndex = startIndex;
        }

        private void combineRead(String otherRead, int otherReadStartIndex) {
            int substringStart = this.startIndex +
                    this.stringBuilder.length() - otherReadStartIndex;
            if (substringStart < otherRead.length()) {
                this.stringBuilder.append(otherRead.substring(substringStart));
            }
        }

        private String build() {
            return this.stringBuilder.toString();
        }
    }






}
