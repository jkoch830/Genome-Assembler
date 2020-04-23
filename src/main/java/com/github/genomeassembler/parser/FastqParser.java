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




    public static List<String> getReads(String path) {
        return new ArrayList<>();
    }
    /**
     * HIV:
     * K = 20 => 17130229
     * K = 30 => 23336643
     * K = 50 => 31383206
     *
     * E.Coli
     *
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

    public static void main(String[] args) {
        findKmers("src/main/resources/genitalium1A.fastq");
        findKmers("src/main/resources/genitalium1B.fastq");
    }

    /*
    Total iterations: 2627780
Total kmers: 2491358
Free memory available: 8132844032
0
1000000
2000000
Total iterations: 2627780
Total kmers: 5388574
Free memory available: 7541765632
     */

}
