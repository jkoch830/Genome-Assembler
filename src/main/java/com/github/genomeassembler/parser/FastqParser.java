package com.github.genomeassembler.parser;


import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

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
        String line;
        for (line = br.readLine(); line != null; line = br.readLine()) {
            if (line.charAt(0) != '@' && nextHasRead) {
                reads.add(line.replaceAll("N", "G"));
                nextHasRead = false;
            } else if (line.charAt(0) == '@') {
                nextHasRead = true;
            }
        }
        return reads;
    }
}
