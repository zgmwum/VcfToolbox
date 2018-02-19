package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Read multi sample vcf file; must be already splitted on multiallelic places
 * 
 * Writes it back, but with additional marks, that enables unique changes
 * identification
 * 
 * @author pstawinski
 * 
 */
public class IdentifySomaticInGroupMpileup {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(IdentifySomaticInGroupMpileup.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    private boolean ungzip = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        IdentifySomaticInGroupMpileup app = new IdentifySomaticInGroupMpileup();
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
        // TODO Auto-generated method stub

        InputStream inputStream;
        if ("-".equals(inputVcfFile)) {
            inputStream = System.in;
        } else {
            inputStream = new FileInputStream(inputVcfFile);
            if (inputVcfFile.endsWith(".gz")) {
                ungzip = true;
            }
        }

        if (ungzip) {
            inputStream = new GZIPInputStream(inputStream);
        }

        VcfStreamReader vcfStreamReader = new VcfStreamReader(new BufferedReader(new InputStreamReader(inputStream)));
        VCFHeader header = vcfStreamReader.getVcfHeader();

        String[] samples = vcfStreamReader.getSamplesLine();

        OutputStream outputStream;
        if ("-".equals(outputVcfFile)) {
            outputStream = System.out;
        } else {
            outputStream = new FileOutputStream(outputVcfFile);
        }
        VcfStreamWriter vcfStreamWriter = new VcfStreamWriter(new BufferedWriter(new OutputStreamWriter(outputStream)));

        header.addMetaDataLine(new VCFInfoHeaderLine("STRONG_PRESENT", 1, VCFHeaderLineType.Integer,
                "count of samples it is present in"));
        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        int counter = 0;

        VcfLine vc;

        while ((vc = vcfStreamReader.next()) != null) {

            Map<String, SampleData> samplesData = new HashMap<>();
            for (String sample : samples) {
                String ads = vc.getSamplesData().get(sample).get("AD");

                if (ads != null && !ads.equals(".")) {
                    String[] adsSplitted = ads.split(",");
                    int adRef = Integer.valueOf(adsSplitted[0]);
                    int adAlt = Integer.valueOf(adsSplitted[1]);
                    SampleData sd = new SampleData();
                    sd.setAltAD(adAlt);
                    sd.setRefAD(adRef);
                    samplesData.put(sample, sd);

                    vc.getSamplesData().get(sample).put("DP", String.valueOf(sd.getDP()));

                }
            }

            int present = samplesData.values().stream().mapToInt(s -> s.presentStrong() ? 1 : 0).sum();
            vc.getInfo().put("STRONG_PRESENT", present);

            vc.getFormat().add("DP");
            vc.setAlt(vc.getAlt().stream().map(String::toUpperCase).collect(Collectors.toList()));
            vc.setRef(vc.getRef().toUpperCase());

            vcfStreamWriter.write(vc);

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfStreamWriter);
        IOUtils.closeQuietly(vcfStreamReader);

    }

    private static class SampleData {
        private String name;
        private int refAD;
        private int altAD;

        public SampleData() {
        }

        public int getRefAD() {
            return refAD;
        }

        public void setRefAD(int refAD) {
            this.refAD = refAD;
        }

        public int getAltAD() {
            return altAD;
        }

        public void setAltAD(int altAD) {
            this.altAD = altAD;
        }

        public int getDP() {
            return refAD + altAD;
        }

        public boolean presentWeakest() {
            return altAD > 0;
        }

        public boolean nodata() {
            return (altAD + refAD) == 0;
        }

        public boolean covered() {
            return (altAD + refAD) >= 10;
        }

        public boolean presentWeak() {
            double altFraction = altAD / (double) (altAD + refAD);
            return altFraction > 0.05 && altAD > 1;
        }

        public boolean presentStrong() {
            double altFraction = altAD / (double) (altAD + refAD);
            return altFraction > 0.1 && altAD > 2;
        }

    }
}
