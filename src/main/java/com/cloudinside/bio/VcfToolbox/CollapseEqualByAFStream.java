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
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;

import htsjdk.variant.vcf.VCFHeader;

/**
 * Adjacent equal variants are connected, AC and AF are summed up
 * 
 * @author pstawinski
 * 
 */
public class CollapseEqualByAFStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(CollapseEqualByAFStream.class);

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

        CollapseEqualByAFStream app = new CollapseEqualByAFStream();
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
                    if (previousLine.getChr().equals(vc.getChr()) && previousLine.getPosition().equals(vc.getPosition())
                            && previousLine.getAlt().equals(vc.getAlt())) {
                        collapse(previousLine, vc, "AF", true);
                        collapse(previousLine, vc, "AC", false);

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

    private void collapse(VcfLine previousLine, VcfLine vc, String field, boolean isDouble) {
        if (previousLine.getInfo().containsKey(field) && vc.getInfo().containsKey(field)) {
            previousLine.getInfo()
                    .put(field, isDouble
                            ? String.format("%.6f",
                                    ((Double.valueOf((String) previousLine.getInfo().get(field))
                                            + Double.valueOf((String) vc.getInfo().get(field)))),
                                    Locale.ENGLISH)
                            : Integer.toString((Integer.valueOf((String) previousLine.getInfo().get(field))
                                    + Integer.valueOf((String) vc.getInfo().get(field)))));
        }
    }
}
