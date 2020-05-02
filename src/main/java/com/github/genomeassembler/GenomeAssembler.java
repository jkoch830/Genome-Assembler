package com.github.genomeassembler;

import com.github.genomeassembler.debruijn.BasicDeBruijnGraph;
import com.github.genomeassembler.debruijn.DeBruijnAnalyzer;
import com.github.genomeassembler.debruijn.DeBruijnGraph;
import com.github.genomeassembler.mapper.BWReadMapper;
import com.github.genomeassembler.mapper.ReadMapper;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class GenomeAssembler {

    // Status messages
    private final static String PROCESS_BEGIN_MSG = "Beginning processing of " +
            "reference genome...";
    private final static String PROCESS_FINISH_MSG = "Finished processing of " +
            "reference genome";
    private final static String MAPPING_MSG = "Mapping reads to reference " +
            "genome with mismatch tolerance: ";
    private final static String MAPPED_MSG = "Number of reads mapped to " +
            "reference genome: ";
    private final static String FORMING_CONTIGS_MSG = "Forming contigs...";
    private final static String CONTIGS_FORMED_MSG = "Number of contigs formed: ";
    private final static String REMAINING_READS_MSG = "Number of remaining " +
            "unmapped reads: ";
    private final static String CONSTRUCTING_GRAPH_MSG = "Constructing de " +
            "Bruijn graph from remaining unmapped reads and assembling contigs...";
    private final static String RESOLVING_MSG = "Resolving and combining contigs";
    private final static String REMOVING_CONTIGS_MSG = "Removing contigs with " +
            "length below threshold: ";
    private final static String REMAINING_CONTIGS_MSG = "Number of remaining " +
            "contigs from graph: ";
    private final static String FINISHED_MSG = "Contig generation " +
            "complete\nWriting to file...";

    private final static int NUM_THREADS = 4;

    private final ReadMapper referenceGenomeReadMapper;
    private final int refGenomeLength;
    private final DeBruijnGraph deBruijnGraph;

    private final List<String> unmappedReads;
    private final Map<String, List<Integer>> mappedReads;

    // Each correspond to a certain mismatch-tolerated mappedReads mapping
    private final List<Map<String, Integer>> mappedContigSets;

    // Super contigs formed by mapped contig sets, then resolved by graph contigs
    private final Map<Integer, String> superContigs;

    // Contigs formed by de Bruijn graph traversal
    private final List<String> graphContigs;

    // Settings
    private AssemblerParameters parameters;

    // Coverage
    private final int coverage;


    public GenomeAssembler(String referenceGenomeSequence, List<String> reads) {
        // Process genome
        this.refGenomeLength = referenceGenomeSequence.length();
        System.out.println("Length of genome: " + referenceGenomeSequence.length());
        System.out.println("Total number of reads: " + reads.size());
        System.out.println(PROCESS_BEGIN_MSG);
        this.referenceGenomeReadMapper = new BWReadMapper(referenceGenomeSequence);
        System.out.println(PROCESS_FINISH_MSG);

        // Initialize read data structures
        this.unmappedReads = reads;
        this.mappedReads = new ConcurrentHashMap<>();

        // Initialize contig data structures
        this.mappedContigSets = new ArrayList<>();
        this.superContigs = new HashMap<>();
        this.graphContigs = new ArrayList<>();

        // Initialize de Bruijn graph
        this.deBruijnGraph = new BasicDeBruijnGraph();

        // Initialize default diagnostic information settings
        this.parameters = new AssemblerParameters.Builder().build();

        // Calculate coverage
        int totalBp = 0;
        for (String read : reads) {
            totalBp += read.length();
        }
        this.coverage = totalBp / referenceGenomeSequence.length();
    }

    /**
     * Main method to assemble the genome after construction.
     * Main process:
     *      - For all x between a given lower bound and higher bound:
     *              - Map all remaining reads that map to reference genome with
     *                  at most x mismatches
     *              - Map all remaining reverse complementary reads that map to
     *                  to reference genome with at most x mismatches
     *              - Form contigs from these mapped reads
     *      - Take the remaining reads and obtain remaining contigs by
     *          constructing a de Bruijn graph. Forms kmers based on configurable
     *          k
     *      - Denote each contig set C-x, where each contig in C-x was formed
     *          by reads that mapped to the reference genome with at most x
     *          mismatches. To form the super contig between multiple
     *          overlapping contigs, take the consensus base at the overlapping
     *          index. To break ties, prioritize the base in C-y, where y is
     *          minimal.
     *      - Resolve gaps between super contigs by mapping graph contigs
     *          onto super contigs. Remove the graph contigs from the set
     *      - From remaining graph contigs, discard contigs that fall below
     *          configurable threshold size
     *      - Output N50, largest contig, super contig(s), graph contigs,
     *          super contig coverage, total contig coverage, and read coverage
     */
    public void assemble() {
        long startTime = System.currentTimeMillis();
        // Use reference genome to form contigs
        int mismatchLowerBound = this.parameters.getMismatchToleranceLowerBound();
        int mismatchUpperBound = this.parameters.getMismatchToleranceHigherBound();
        int toleranceStep = this.parameters.getMismatchToleranceStep();
        int mapped, contigs;
        for (int numMismatches = mismatchLowerBound;
             numMismatches < mismatchUpperBound; numMismatches += toleranceStep) {
            System.out.println(MAPPING_MSG + numMismatches);
            mapped = this.mapReads(numMismatches); // Tries mapping reads
            System.out.println(MAPPED_MSG + mapped);
            System.out.println(FORMING_CONTIGS_MSG);
            contigs = this.formContigs(); // Forms contigs from reads
            System.out.println(CONTIGS_FORMED_MSG + contigs);
        }
        System.out.println(REMAINING_READS_MSG + this.unmappedReads.size());
        System.out.println(CONSTRUCTING_GRAPH_MSG);

        // Form contigs out of remaining reads using a de Bruijn graph
        System.out.println();
        constructDeBruijnGraph();
        this.graphContigs.addAll(DeBruijnAnalyzer.contigGeneration(this.deBruijnGraph));
        System.out.println(CONTIGS_FORMED_MSG + this.graphContigs.size());
        System.out.println(REMOVING_CONTIGS_MSG + this.parameters.getMinContigOutputLength());
        this.removeSmallGraphContigs();
        System.out.println(REMAINING_CONTIGS_MSG + this.graphContigs.size());

        // Map graph contigs to reference and form contigs
        this.unmappedReads.addAll(this.graphContigs);
        this.mapReads(0);
        this.formContigs();

        // Combine/resolve contigs
        System.out.println(RESOLVING_MSG);
        this.resolveContigs();

        // Use super contigs to write out assembled genome sequence with gaps
        String assembledGenome = this.formSequence();

        // Write assembly results in file
        System.out.println(FINISHED_MSG);
        long endTime = System.currentTimeMillis();
        this.writeResults(assembledGenome, endTime - startTime);

    }

    /* *************************************************** */
    /* ************** BEGIN WRITING METHODS ************** */
    /* *************************************************** */
    /**
     * Outputs N50, largest contig length, coverage, and contigs
     */
    private void writeResults(String assembledGenome, long totalTime) {
        int N50 = this.calculateN50();
        int longestContigLength = this.getLongestContigSize();
        double allContigCoverage = this.calculateAllContigCoverage();
        double superContigCoverage = this.calculateSuperContigCoverage();
        int numGaps = countChar(assembledGenome, 'N');
        List<Integer> sortedIndices = getSortedIndices(this.superContigs);
        try {
            String path = "src/main/resources/results.txt";
            BufferedWriter writer = new BufferedWriter(new FileWriter(path));
            writer.write("PARAMETERS:\n");
            writer.write("k: " + this.parameters.getKmerLength() + "\n");
            writer.write("Mismatch tolerance lower bound: " +
                    this.parameters.getMismatchToleranceLowerBound() + "\n");
            writer.write("Mismatch tolerance higher bound: " +
                    this.parameters.getMismatchToleranceHigherBound() + "\n");
            writer.write("Mismatch tolerance step: " +
                    this.parameters.getMismatchToleranceStep() + "\n");
            writer.write("\n");
            writer.write("ASSEMBLY DETAILS:\n");
            writer.write("Total time: " + (totalTime / 60000) + " minutes\n");
            writer.write("N50: " + N50 + "\n");
            writer.write("Longest Contig Length: " + longestContigLength + "bp\n");
            writer.write("Coverage of reads: " + this.coverage + "x\n");
            writer.write("Coverage of super contigs: " + superContigCoverage + "%\n");
            writer.write("Coverage of contigs: " + allContigCoverage + "%\n");
            writer.write("Number of unknown bases: " + numGaps + "\n");
            writer.write("\n");
            writer.write("======================Super Contigs======================\n");
            for (int index : sortedIndices) {
                String contig = this.superContigs.get(index);
                int size = contig.length();
                writer.write("Super Contig index: " + index + "\n");
                writer.write("Super Contig length: " + size + "\n");
                writer.write(contig + "\n");
            }
            writer.write("\n");
            writer.write("======================Graph Contigs======================\n");
            for (String contig : this.graphContigs) {
                int size = contig.length();
                writer.write("Contig length: " + size + "\n");
                writer.write(contig + "\n");
            }
            writer.write("\n");
            writer.write("======================Assembled Genome======================\n");
            writer.write(assembledGenome + "\n");
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error writing to file");
            System.out.println("N50: " + N50);
            System.out.println("Longest contig: " + longestContigLength);
            System.out.println("Coverage: " + coverage);
            System.out.println("Super Contig Coverage: " + superContigCoverage);
            System.out.println("All Contig Coverage: " + allContigCoverage);
        }
    }

    /**
     * Uses the super contigs to form an assembled sequence. Unknown regions
     * are written by N's
     * @return The assembled sequence
     */
    private String formSequence() {
        StringBuilder assembledSequence = new StringBuilder();
        List<Integer> sortedIndices = getSortedIndices(this.superContigs);
        int index = 0;
        while (index < this.refGenomeLength) {
            if (!sortedIndices.contains(index)) { // Gap
                assembledSequence.append('N');
                index++;
            } else {
                String superContig = this.superContigs.get(index);
                assembledSequence.append(superContig);
                index += superContig.length();
            }
        }
        return assembledSequence.toString();
    }

    private int getLongestContigSize() {
        int longestContigLength = 0;
        for (String superContig : this.superContigs.values()) {
            if (superContig.length() > longestContigLength) {
                longestContigLength = superContig.length();
            }
        }
        return longestContigLength;
    }

    private double calculateSuperContigCoverage() {
        int totalLength = 0;
        for (String contig : this.superContigs.values()) {
            totalLength += contig.length();
        }
        return 100 * totalLength / ((double)this.refGenomeLength);
    }

    private double calculateAllContigCoverage() {
        int totalLength = 0;
        for (String superContig : this.superContigs.values()) {
            totalLength += superContig.length();
        }
        for (String contig : this.graphContigs) {
            totalLength += contig.length();
        }
        return 100 * totalLength / ((double)this.refGenomeLength);
    }

    private int calculateN50() {
        List<Integer> sizes = new ArrayList<>();
        for (String superContig : this.superContigs.values()) {
            sizes.add(superContig.length());
        }
        Collections.sort(sizes);
        Collections.reverse(sizes);
        int total = 0;
        for (int size : sizes) {
            total += size;
            if (total >= (this.refGenomeLength/ 2)) {
                System.out.println("Good N50 calculation: " + size);
                return size;
            }
        }
        System.out.println("Bad n50 Calculation");
        System.out.println("Total only reached: " + total);
        return total;
    }

    /* ************************************************* */
    /* ************** END WRITING METHODS ************** */
    /* ************************************************* */



    /* *************************************************** */
    /* ************** BEGIN ASSEMBLER STEPS ************** */
    /* *************************************************** */
    /**
     * Resolves and combines contigs from the mapped contig sets.
     * It forms consensus super contigs by taking the consensus base from each
     * mapped contig.
     * This places the super contigs and clears the mapped contig sets
     */
    private void resolveContigs() {
        // Flips all contig sets
        List<Map<Integer, String>> flippedContigSets = new ArrayList<>();
        for (Map<String, Integer> contigSet : this.mappedContigSets) {
            Map<Integer, String> flippedContigSet = new HashMap<>();
            for (Map.Entry<String, Integer> mapping : contigSet.entrySet()) {
                flippedContigSet.put(mapping.getValue(), mapping.getKey());
            }
            flippedContigSets.add(flippedContigSet);
            contigSet.clear();
        }
        // Begin forming consensus contig between mapped contigs
        StringBuilder superContig = new StringBuilder();
        for (int i = 0; i < this.refGenomeLength; i++) {
            char consensusBase = getConsensusBase(i, flippedContigSets);
            if (consensusBase == '0') { // No contig covers i
                if (superContig.length() != 0) { // don't save empty string
                    int startingIndex =  i - superContig.length();
                    this.superContigs.put(startingIndex, superContig.toString());
                    superContig = new StringBuilder();
                }
            } else {
                superContig.append(consensusBase);
            }
        }
        // Checks for last super contig
        int finalSuperContigLength = superContig.length();
        if (finalSuperContigLength != 0) {
            int startingIndex = this.refGenomeLength - finalSuperContigLength;
            this.superContigs.put(startingIndex, superContig.toString());
        }

        this.mappedContigSets.clear(); // No longer need contigs
    }


    /**
     * Removes graph contigs that fall below threshold size
     */
    private void removeSmallGraphContigs() {
        int maxLength = 0;

        int threshold = this.parameters.getMinContigOutputLength();
        for (int i = this.graphContigs.size() - 1; i >= 0; i--) {
            int length = this.graphContigs.get(i).length();
            if (length < threshold) {
                this.graphContigs.remove(i);
            }
            if (length > maxLength) { maxLength = length; }
        }
        System.out.println("BIGGEST CONTIG IN GRAPH: " + maxLength);
    }



    /**
     * Creates a de Bruijn graph from the remaining unmapped reads
     */
    private void constructDeBruijnGraph() {
        int k = this.parameters.getKmerLength();
        for (String read : this.unmappedReads) {
            for (int i = 0; i < read.length() - k + 1; i++) {
                this.deBruijnGraph.addKmer(read.substring(i, i + k));
            }
        }
        this.unmappedReads.clear();
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
        List<Integer> mappedReadsIndices = Collections.synchronizedList(new ArrayList<>());
        int n = this.unmappedReads.size();
        int percentUpdateIncrement = (n < 20) ? 1 : n / 20;

        System.out.println("Mapping " + n + " reads to reference genome...");

        AtomicInteger countRead = new AtomicInteger();
        Thread[] threads = new Thread[NUM_THREADS];
        for (int i = 0; i < threads.length; i++) {
            int startIndex = i * (n / threads.length);
            int endIndex;
            if (i == threads.length - 1) { endIndex = n; }
            else { endIndex = startIndex + (n / threads.length); }
            threads[i] = new Thread(() -> {
                for (int j = startIndex; j < endIndex; j++) {
                    countRead.getAndIncrement();
                    if (countRead.get() % percentUpdateIncrement == 0) {
                        int percentComplete = 5 * (countRead.get() / percentUpdateIncrement);
                        System.out.println(percentComplete + "% reads mapped");
                    }
                    String read = this.unmappedReads.get(j);
                    List<Integer> startingPositions =
                            this.referenceGenomeReadMapper.mapRead(read, mismatches);
                    if (startingPositions.isEmpty()) { // Try mapping complement
                        read = reverseComplement(read);
                        startingPositions = this.referenceGenomeReadMapper.mapRead(
                                read, mismatches);

                    }
                    // Add if and only if it was mapped to at least one position
                    if (!startingPositions.isEmpty()) {
                        this.mappedReads.put(read, startingPositions);
                        mappedReadsIndices.add(j);
                    }
                }
            });
        }
        for (Thread thread : threads) { thread.start(); }
        for (Thread thread : threads) {
            try { thread.join(); }
            catch (InterruptedException e) { e.printStackTrace(); }
        }
        Collections.sort(mappedReadsIndices);
        Collections.reverse(mappedReadsIndices);
        for (int index : mappedReadsIndices) { // Safely removes all reads that were mapped
            this.unmappedReads.remove(index);
        }
        return mappedReadsIndices.size();
    }


    /**
     * Forms contigs from mapped reads that form a contiguous sequence.
     * The mapping containing the reads will be cleared because the resulting
     * contig map will contain reads that don't form a contig. The new contig
     * mapping will be added to the assembler's mapped contig set list
     * @return The number of contigs formed
     */
    public int formContigs() {
        if (this.mappedReads.isEmpty()) {
            return 0;
        }
        Map<Integer, String> mappedIndices = removeSmallStrings(flipMapping(this.mappedReads));
        this.mappedReads.clear(); // No need to store reads any longer

        Map<String, Integer> contigs = new HashMap<>(); // New contig set

        // Traverse mapped contig positions in order
        Set<Integer> tempKeySet = mappedIndices.keySet();
        int start = Collections.min(tempKeySet);
        int end = Collections.max(tempKeySet);
        tempKeySet = null;

        int requiredOverlap = this.parameters.getRequiredContigOverlap(), readLength;
        ContigBuilder newContig;
        for (int i = start; i <= end; i++) { // Loop through mapped indices
            if (mappedIndices.containsKey(i)) { // Found read at index
                String read = mappedIndices.get(i);
                readLength = read.length();
                newContig = new ContigBuilder(read, i);
                // Tries finding other strings it overlaps with
                int j = i + 1;
                int validSearchBound = i + readLength - requiredOverlap, newBound;
                while (j <= validSearchBound) {
                    if (mappedIndices.containsKey(j)) { // overlaps with other read
                        String otherRead = mappedIndices.get(j);
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
        this.mappedContigSets.add(contigs); // Add new contig set
        return contigs.size();
    }

    /* ************************************************* */
    /* ************** END ASSEMBLER STEPS ************** */
    /* ************************************************* */

    /**
     * Sets the assembler's settings
     * @param parameters An AssemblerParameters object containing assembly
     *                   settings
     */
    public void setAssemblerParameters(AssemblerParameters parameters) {
        this.parameters = parameters;
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
     * @return A map containing each mapped read and their respective starting
     *      positions in the reference genome
     */
    public Map<String, List<Integer>> getMappedReads() {
        return this.mappedReads;
    }


    /**
     * Retrieve's the assembler's already formed contig sets
     * @return A mapping between the contigs and their starting position in the
     *         genome
     */
    public List<Map<String, Integer>> getMappedContigSets() {
        return this.mappedContigSets;
    }



    /** Static helper methods and classes */
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

    private static List<Integer> getSortedIndices(Map<Integer, String> map) {
        List<Integer> indices = new ArrayList<>(map.keySet());
        Collections.sort(indices);
        return indices;
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

    private static int countChar(String s, char target) {
        int total = 0;
        for (char c : s.toCharArray()) {
            if (c == target) {
                total++;
            }
        }
        return total;
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
