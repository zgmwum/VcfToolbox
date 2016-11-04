package com.cloudinside.bio.VcfToolbox.singleshot;

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
import java.util.Collections;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class SplitMultiallelic1000genomesVcf {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(SplitMultiallelic1000genomesVcf.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = false)
    private String inputVcfFile = "/wum/pio/tmp/1000g/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz";
    @Parameter(names = "--output", description = "Vcf output", required = false)
    private String outputVcfFile = "/tmp/1000g.vcf";

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    private DecimalFormat df = new DecimalFormat("0.00000");

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        SplitMultiallelic1000genomesVcf app = new SplitMultiallelic1000genomesVcf();
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

        System.err.println("Input: " + inputVcfFile);
        System.err.println("Output: " + outputVcfFile);

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

        int counter = 0;

        VcfLine vc;
        while ((vc = vcfStreamReader.next()) != null) {
            if (vc.getAlt().size() == 1) {
                vcfStreamWriter.write(vc);
            } else {
                for (int i = 0; i < vc.getAlt().size(); i++) {
                    VcfLine part = new VcfLine(vc);
                    part.setAlt(Collections.singletonList(vc.getAlt().get(i)));
                    for (Entry<String, Object> info : vc.getInfo().entrySet()) {
                        String val = (String) info.getValue();
                        if (val.contains(",")) {
                            if (val.split(",").length == vc.getAlt().size()) {
                                part.getInfo().put(info.getKey(), val.split(",")[i]);
                            }
                        }
                    }
                    vcfStreamWriter.write(part);
                }
            }

            if (++counter % 10000 == 0) {
                // System.err.println("Processed " + cou. nter);
            }
        }

        IOUtils.closeQuietly(vcfStreamWriter);
        IOUtils.closeQuietly(vcfStreamReader);

    }
}
