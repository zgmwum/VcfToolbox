package com.cloudinside.bio.VcfToolbox.helpers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

/**
 * Can collapse rows with different EFF_0__AA and EFF_0__TRID only to single row
 * 
 * @author pstawinski
 * 
 */
public class CollapseEffTranscriptRows {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(CollapseEffTranscriptRows.class);

    @Parameter(names = "--input", description = "Tsv input file ( - for stdin)", required = true)
    private String inputFile;
    @Parameter(names = "--output", description = "Tsv output ( - for stdout)", required = true)
    private String outputFile;

    // @Parameter(names = "--joinEffAndImpact", description =
    // "Join effect and impact fields", required = true)
    // private String outputFile;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        CollapseEffTranscriptRows app = new CollapseEffTranscriptRows();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private void go() throws IOException {
        BufferedReader br;
        if (inputFile.equals("-")) {
            br = new BufferedReader(new InputStreamReader(System.in));
        } else {
            br = new BufferedReader(new FileReader(inputFile));
        }

        BufferedWriter bw;
        if (outputFile.equals("-")) {
            bw = new BufferedWriter(new OutputStreamWriter(System.out));
        } else {
            bw = new BufferedWriter(new FileWriter(outputFile));
        }

        Splitter splitter = Splitter.on('\t');
        Joiner joiner = Joiner.on('\t');
        Joiner joinerSemicolon = Joiner.on(';');
        String line;
        List<String> header = new ArrayList<>(splitter.splitToList(br.readLine()));
        int aaChangeIndex = header.indexOf("EFF_0__AA");
        int tridIndex = header.indexOf("EFF_0__TRID");

        header.remove("EFF_0__AA");
        header.remove("EFF_0__TRID");
        int newIndex = Math.min(aaChangeIndex, tridIndex);
        header.add(newIndex, "EFF_AA");

        bw.write(joiner.join(header) + "\n");

        List<String> previousLine = null;
        List<String> aaSeen = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            List<String> lineSplitted = new ArrayList<>(splitter.splitToList(line));
            String currentAA = lineSplitted.get(aaChangeIndex);
            String currentTrid = lineSplitted.get(tridIndex);
            String currentAATrid = currentTrid + ":" + currentAA;
            if (currentAATrid.equals(":"))
                currentAATrid = "";
            String removed = lineSplitted.remove(Math.min(aaChangeIndex, tridIndex));
            removed = lineSplitted.remove(Math.max(aaChangeIndex, tridIndex) - 1);
            if (previousLine == null) {
                previousLine = lineSplitted;
                aaSeen.add(currentAATrid);
            } else if (equal(previousLine, lineSplitted)) {
                aaSeen.add(currentAATrid);
            } else {
                previousLine.add(newIndex, joinerSemicolon.join(aaSeen));
                bw.write(joiner.join(previousLine) + "\n");
                aaSeen.clear();
                previousLine = lineSplitted;
                aaSeen.add(currentAATrid);
            }
            // System.err.println(line);

        }
        previousLine.add(newIndex, joinerSemicolon.join(aaSeen));
        bw.write(joiner.join(previousLine) + "\n");

        br.close();
        bw.close();
    }

    private boolean equal(List<String> previousLine, List<String> lineSplitted) {
        for (int i = 0; i < previousLine.size(); i++) {
            if (!previousLine.get(i).equals(lineSplitted.get(i)))
                return false;
        }
        return true;
    }
}
