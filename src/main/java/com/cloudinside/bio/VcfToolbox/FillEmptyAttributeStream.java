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
import java.util.List;
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;

import htsjdk.variant.vcf.VCFHeader;

/**
 * Finds variants in genes that are possibly damaging and appears twice (i.e.
 * possibly two gene copies are destroyed).
 * 
 * @author pstawinski
 * 
 */
public class FillEmptyAttributeStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(FillEmptyAttributeStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix). Cann vcf file with snpEff annotations", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    @Parameter(names = "--attribute", description = "Attribute name", required = true)
    private List<String> attributes;
    @Parameter(names = "--value", description = "Attribute value", required = true)
    private List<String> values;
    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        FillEmptyAttributeStream app = new FillEmptyAttributeStream();
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

        if (ungzip) {
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
            VcfLine vc;
            while ((vc = it.next()) != null) {
                int index = 0;
                for (String attribute : attributes) {
                    if (vc.getInfo() != null) {
                        Object attributeValue = vc.getInfo().get(attribute);
                        if (attributeValue == null) {
                            vc.getInfo().put(attribute, values.get(index));
                        }
                        index++;
                    } else {
                        log.error("I've received empty info for:" + vc);
                    }
                }

                vcfStreamWriter.write(vc);

            }
        }

        vcfStreamReader.close();
        vcfStreamWriter.close();

    }
}
