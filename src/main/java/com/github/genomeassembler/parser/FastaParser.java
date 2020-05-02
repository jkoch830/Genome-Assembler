package com.github.genomeassembler.parser;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Class containing methods used to parse genome files
 */
public class FastaParser {
    public static String parseGffFile(String path) {
        int iteration = 0;
        boolean atSequence = false;
        StringBuilder genome = new StringBuilder();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8)) {
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                if (iteration % 1000000 == 0) {
                    System.out.println(iteration);
                }
                iteration++;
                if (atSequence) {
                    genome.append(line.replaceAll("\\s", ""));
                } else if (line.substring(0, 7).equals("##FASTA")){
                    // Waits until it reaches fasta sequence
                    br.readLine(); // Skip next line
                    atSequence = true;
                }
            }
        } catch (IOException e) {
            System.out.println("Error with reading file");
            e.printStackTrace();
        }
        return genome.toString();
    }
}
