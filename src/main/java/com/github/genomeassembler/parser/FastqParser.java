package com.github.genomeassembler.parser;


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
}
