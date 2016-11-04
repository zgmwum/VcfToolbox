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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.broadinstitute.variant.vcf.VCFHeader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;
import com.google.common.collect.ImmutableSet;

/**
 * Adjacent equal variants are connected, AC and AF are summed up; can be used
 * for multisample vcf; cannot be used for multi-observation vcf
 * 
 * @author pstawinski
 * 
 */
public class CollapseEqualMultisampleStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(CollapseEqualMultisampleStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix). Cann vcf file with snpEff annotations", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    private NumberFormat formatter = new DecimalFormat("#0.00");

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        CollapseEqualMultisampleStream app = new CollapseEqualMultisampleStream();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private final static int TRANSCRIPT_NAME_POSITION = 10;
    private final static Set<String> differentFieldsThatCanBeOmmited = ImmutableSet.copyOf(new String[] { "AF", "AC" });

    private void go() throws IOException {
        // TODO Auto-generated method stub
        InputStream inputStream;
        if ("-".equals(inputVcfFile)) {
            inputStream = System.in;
        } else {
            inputStream = new FileInputStream(inputVcfFile);
        }

        if (ungzip || inputVcfFile.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
        }

        VcfStreamReader vcfStreamReader = new VcfStreamReader(new BufferedReader(new InputStreamReader(inputStream)));
        VCFHeader header = vcfStreamReader.getVcfHeader();

        OutputStream outputStream;
        if ("-".equals(outputVcfFile)) {
            outputStream = System.out;
        } else {
            outputStream = new FileOutputStream(outputVcfFile);
        }
        VcfStreamWriter vcfStreamWriter = new VcfStreamWriter(new BufferedWriter(new OutputStreamWriter(outputStream)));

        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        VcfStreamReader it = vcfStreamReader;

        {
            VcfLine previousLine = null;
            VcfLine vc;

            while ((vc = it.next()) != null) {

                int index = 0;
                if (previousLine == null) {
                    previousLine = vc;
                } else {
                    if (previousLine.getChr().equals(vc.getChr())
                            && previousLine.getPosition().equals(vc.getPosition())
                            && previousLine.getAlt().equals(vc.getAlt())) {
                        collapse(previousLine, vc, "AF", true);
                        collapse(previousLine, vc, "AC", false);
                        deleteDifferentInfoFields(previousLine, vc);
                        mergeGenotypes(previousLine, vc);

                    } else {
                        vcfStreamWriter.write(previousLine);
                        previousLine = vc;
                    }
                }

            }
            if (previousLine != null)
                vcfStreamWriter.write(previousLine);
        }

        vcfStreamReader.close();
        vcfStreamWriter.close();

    }

    private void mergeGenotypes(VcfLine previousLine, VcfLine vc) {
        // I assume identical samples order and GT being in first FORMAT
        // position

        try {
            previousLine.setFormat(Collections.singletonList("GT"));

            LinkedHashMap<String, List<String>> samplesData = previousLine.getSamplesData();
            for (Entry<String, List<String>> entry : samplesData.entrySet()) {
                String previousLineGT = entry.getValue().get(0);
                String vcGT = vc.getSamplesData().get(entry.getKey()).get(0);

                // normalize to make comparisions easier
                previousLineGT = previousLineGT.replace("|", "/").replace("./.", ".");
                vcGT = vcGT.replace("|", "/").replace("./.", ".");

                String gt = null;
                if (previousLineGT.equals(vcGT)) {
                    gt = previousLineGT;
                } else {
                    if ("1/1".equals(previousLineGT) || "1".equals(previousLineGT)) {
                        gt = previousLineGT;
                    } else if ("1/1".equals(vcGT) || "1".equals(vcGT)) {
                        gt = vcGT;
                    } else if ("0/1".equals(previousLineGT)) {
                        gt = previousLineGT;
                    } else if ("0/1".equals(vcGT)) {
                        gt = vcGT;
                    } else if ("0/0".equals(previousLineGT) || "0".equals(previousLineGT)) {
                        gt = previousLineGT;
                    } else if ("0/0".equals(vcGT) || "0".equals(vcGT)) {
                        gt = vcGT;
                    } else if ("./1".equals(previousLineGT)) {
                        gt = previousLineGT;
                    } else if ("./1".equals(vcGT)) {
                        gt = vcGT;
                    } else {
                        System.err.println("in one place we have: " + vcGT + " and " + previousLineGT + " at "
                                + vc.getChr() + ":" + vc.getPosition());
                    }
                }
                entry.setValue(Collections.singletonList(gt));
            }
        } catch (Exception e) {
            System.err.println("Exception at " + vc.getChr() + ":" + vc.getPosition());
            System.err.println(e.getMessage());
            System.err.println(vc);
            e.printStackTrace();
        }
    }

    private void deleteDifferentInfoFields(VcfLine previousLine, VcfLine vc) {
        Iterator<Entry<String, Object>> it = previousLine.getInfo().entrySet().iterator();
        while (it.hasNext()) {
            Entry<String, Object> entry = it.next();
            if (!differentFieldsThatCanBeOmmited.contains(entry.getKey())) {
                if (entry.getValue().equals(vc.getInfo().get(entry.getKey()))) {
                    // leave
                } else {
                    it.remove();
                }
            }
        }
    }

    private void collapse(VcfLine previousLine, VcfLine vc, String field, boolean isDouble) {
        if (previousLine.getInfo().containsKey(field) && vc.getInfo().containsKey(field)) {
            previousLine.getInfo().put(
                    field,
                    isDouble ? String.format("%.6f",
                            ((Double.valueOf((String) previousLine.getInfo().get(field)) + Double.valueOf((String) vc
                                    .getInfo().get(field)))), Locale.ENGLISH) : Integer.toString((Integer
                            .valueOf((String) previousLine.getInfo().get(field)) + Integer.valueOf((String) vc
                            .getInfo().get(field)))));
        }
    }
}
