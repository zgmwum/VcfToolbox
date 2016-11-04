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
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VariantHashCounter;
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
public class AddZgmIdStream {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(AddZgmIdStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        AddZgmIdStream app = new AddZgmIdStream();
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

        header.addMetaDataLine(new VCFInfoHeaderLine(ZGM_HASH_INFO_NAME, 1, VCFHeaderLineType.String,
                "unique identifier of this variant"));
        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        int counter = 0;

        VcfLine vc;
        while ((vc = vcfStreamReader.next()) != null) {

            vc.getInfo().put(ZGM_HASH_INFO_NAME,
                    VariantHashCounter.hash(vc.getChr(), vc.getPosition(), null, vc.getRef(), vc.getAlt().get(0)));

            vcfStreamWriter.write(vc);

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfStreamWriter);
        IOUtils.closeQuietly(vcfStreamReader);

    }
}
