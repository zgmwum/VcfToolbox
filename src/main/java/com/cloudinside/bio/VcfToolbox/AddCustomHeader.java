package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

/**
 * For adding custom header, f.ex. ##zgm_phenotype=HI
 * 
 * @author pstawinski
 * 
 */
public class AddCustomHeader {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(AddCustomHeader.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--override", description = "Remove all headers with this name before adding", required = false)
    private boolean override = false;
    @Parameter(names = "--key", description = "key", required = true)
    private String key;
    @Parameter(names = "--value", description = "value", required = true)
    private String value;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        AddCustomHeader app = new AddCustomHeader();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private final static String NEWLINE = "\n";

    private void go() {
        // TODO Auto-generated method stub

        BufferedReader reader = null;
        BufferedWriter writer = null;
        try {
            reader = new BufferedReader(new FileReader(new File(inputVcfFile)));
            writer = new BufferedWriter(new FileWriter(new File(outputVcfFile)));

            String line = null;

            // copy and modify header
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) {
                    if (!override || !line.startsWith("##" + key + "=")) {
                        IOUtils.write(line + NEWLINE, writer);
                    } else {
                        // do not copy
                    }
                } else if (line.startsWith("#CHROM")) {
                    break;
                }
            }

            // add new header
            IOUtils.write("##" + key + "=" + value + NEWLINE, writer);

            // line contains '#CHROM ...' here
            IOUtils.write(line + NEWLINE, writer);

            // copy rest of file
            IOUtils.copy(reader, writer);

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        IOUtils.closeQuietly(reader);
        IOUtils.closeQuietly(writer);

    }
}
