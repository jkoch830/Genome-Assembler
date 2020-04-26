package com.github.genomeassembler.parser;

import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TLongHashSet;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

public class FastqParser {


    /**
     * Parses a fastq file for all reads
     * @param path The file path
     * @return A list of all reads
     */
    public static List<String> getReads(String path) throws IOException {
        List<String> reads = new ArrayList<>();
        BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8);
        boolean nextHasRead = false;
        String line, read;
        int numN = 0;
        for (line = br.readLine(); line != null; line = br.readLine()) {
            if (line.charAt(0) != '@' && nextHasRead) {
                read = line.replaceAll("\\s", "");
                if (line.contains("N")) {
                    numN++;
                }
                reads.add(line.replaceAll("N", "G"));

                nextHasRead = false;
            } else if (line.charAt(0) == '@') {
                nextHasRead = true;
            }
        }
        System.out.println(numN + " reads had an N");
        return reads;
    }



    /**
     * @param path The file path to the fasta file
     */
    public static void findKmers(String path) {
        THashSet<String> kmers = new THashSet<>();
        int distinctCount = 0, length;
        int iteration = 0;
        int totalKmers = 0;
        int totalReads = 0;
        try (BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8)) {
            boolean nextHasRead = false;
            for (String line = null; (line = br.readLine()) != null;) {
                if (iteration % 1000000 == 0) {
                    System.out.println(iteration);
                }
                if (line.charAt(0) != '@' && nextHasRead) {
                    length = line.length();
                    totalReads++;
                    for (int i = 0; i < length - 50 + 1; i++) {
                        String kmer = line.substring(i, i + 50);
                        if (!kmers.contains(kmer)) {
                            kmers.add(kmer);
                            distinctCount++;
                        }
                        totalKmers++;
                    }
                    nextHasRead = false;
                }
                else if (line.charAt(0) == '@') {
                    nextHasRead = true;
                }
                iteration++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        long usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        long freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
        System.out.println("Total kmers: " + totalKmers);
        System.out.println("Total different kmers: " + distinctCount);
        System.out.println("Total reads: " + totalReads);
        System.out.println("Free memory available: " + freeMemory);


    }

}
