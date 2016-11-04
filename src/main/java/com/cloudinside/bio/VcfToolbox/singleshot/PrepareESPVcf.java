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
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

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
public class PrepareESPVcf {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PrepareESPVcf.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = false)
    private String inputVcfFile = "/tmp/esp/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.modified.vcf";
    @Parameter(names = "--output", description = "Vcf output", required = false)
    private String outputVcfFile = "/tmp/esp/x.vcf";

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    private DecimalFormat df = new DecimalFormat("0.00000");

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        PrepareESPVcf app = new PrepareESPVcf();
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

        header.addMetaDataLine(new VCFInfoHeaderLine("EUR_AF", 1, VCFHeaderLineType.Float,
                "this variant frequency in european population (not MAF)"));
        header.addMetaDataLine(new VCFInfoHeaderLine("AFR_AF", 1, VCFHeaderLineType.Float,
                "this variant frequency in african population (not MAF)"));
        header.addMetaDataLine(new VCFInfoHeaderLine("TOT_AF", 1, VCFHeaderLineType.Float,
                "this variant frequency in total population (not MAF)"));

        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        int counter = 0;

        VcfLine vc;
        while ((vc = vcfStreamReader.next()) != null) {
            String ea_ac = (String) vc.getInfo().get("EA_AC");
            int allele_1_EurAc = Integer.valueOf(ea_ac.split(",")[0]).intValue();
            int allele_2_EurAc = Integer.valueOf(ea_ac.split(",")[1]).intValue();
            String aa_ac = (String) vc.getInfo().get("AA_AC");
            int allele_1_AfrAc = Integer.valueOf(aa_ac.split(",")[0]).intValue();
            int allele_2_AfrAc = Integer.valueOf(aa_ac.split(",")[1]).intValue();

            vc.getInfo().put("EUR_AF", df.format((double) allele_1_EurAc / (double) (allele_1_EurAc + allele_2_EurAc)));
            vc.getInfo().put("AFR_AF", df.format((double) allele_1_AfrAc / (double) (allele_1_AfrAc + allele_2_AfrAc)));
            vc.getInfo().put(
                    "TOT_AF",
                    df.format((double) (allele_1_EurAc + allele_1_AfrAc)
                            / (double) (allele_1_AfrAc + allele_2_AfrAc + allele_1_EurAc + allele_2_EurAc)));

            vcfStreamWriter.write(vc);

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfStreamWriter);
        IOUtils.closeQuietly(vcfStreamReader);

    }
}
