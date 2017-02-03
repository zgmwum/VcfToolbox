/**
 * 
 */
package com.cloudinside.bio.VcfToolbox.singleshot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * In hg-summary.vcf
 * 
 * @author pstawinski
 * 
 */
public class FindRecessiveChanges {

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    private final static int GENE_INDEX = 3;
    private final static int SAMPLES_INDEX = 40;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        FindRecessiveChanges app = new FindRecessiveChanges();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private final Splitter commaSplitter = Splitter.on(',');
    private final Splitter pipeSplitter = Splitter.on('|');
    private final Joiner commaJoiner = Joiner.on(',');

    public void go() {
        try {
            Multiset<String> genesToSample = HashMultiset.create();

            InputStream inputStream;

            inputStream = new FileInputStream(inputVcfFile);

            if (inputVcfFile.endsWith(".gz")) {
                inputStream = new GZIPInputStream(inputStream);
            }

            VcfStreamReader vcfStreamReader = new VcfStreamReader(
                    new BufferedReader(new InputStreamReader(inputStream)));
            VCFHeader header = vcfStreamReader.getVcfHeader();

            OutputStream outputStream;
            if ("-".equals(outputVcfFile)) {
                outputStream = System.out;
            } else {
                outputStream = new FileOutputStream(outputVcfFile);
            }

            header.addMetaDataLine(new VCFInfoHeaderLine("ZGM_Samples_Recessive", 1, VCFHeaderLineType.String,
                    "samples having two variants in one gene"));

            VcfStreamWriter vcfStreamWriter = new VcfStreamWriter(
                    new BufferedWriter(new OutputStreamWriter(outputStream)));

            vcfStreamWriter.setVcfHeader(header);
            vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

            VcfLine line;

            while ((line = vcfStreamReader.next()) != null) {

                Set<String> genes = extractGeneNames(line);
                Set<String> samples = extractSamples(line);

                for (String sample : samples) {
                    for (String gene : genes) {
                        sample = sample.replace("*", "");
                        genesToSample.add(gene + ">>" + sample);
                    }
                }
            }
            vcfStreamReader.close();

            // Open again
            inputStream = new FileInputStream(inputVcfFile);

            if (inputVcfFile.endsWith(".gz")) {
                inputStream = new GZIPInputStream(inputStream);
            }

            vcfStreamReader = new VcfStreamReader(new BufferedReader(new InputStreamReader(inputStream)));
            // second pass

            while ((line = vcfStreamReader.next()) != null) {
                List<String> samplesPassed = new ArrayList<>();

                Set<String> genes = extractGeneNames(line);
                Set<String> samples = extractSamples(line);

                for (String sample : samples) {
                    for (String gene : genes) {
                        sample = sample.replace("*", "");
                        int count = genesToSample.count(gene + ">>" + sample);
                        if (count > 1) {
                            samplesPassed.add(sample);

                        }
                    }
                }

                if (!samplesPassed.isEmpty()) {
                    Collections.sort(samplesPassed);
                    line.getInfo().put("ZGM_Samples_Recessive", commaJoiner.join(samplesPassed));
                    vcfStreamWriter.write(line);
                }
            }
            vcfStreamWriter.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    private Set<String> extractSamples(VcfLine line) {
        String eff = (String) line.getInfo().get("ZGM_Samples");
        List<String> samples = commaSplitter.splitToList(eff);
        return new HashSet<String>(samples);
    }

    private Set<String> extractGeneNames(VcfLine line) {
        Set<String> geneNames = new HashSet<>();
        String eff = (String) line.getInfo().get("EFF");
        List<String> effects = commaSplitter.splitToList(eff);

        for (String effect : effects) {
            effect = StringUtils.removeEnd(StringUtils.removeStart(effect, "("), ")");
            geneNames.add(pipeSplitter.splitToList(effect).get(5));
        }
        return geneNames;
    }
}
