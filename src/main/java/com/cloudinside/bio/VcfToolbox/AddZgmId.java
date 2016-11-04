package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.Locale;

import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VariantHashCounter;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class AddZgmId {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(AddZgmId.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header", required = false)
    private boolean ignoreMissingHeader = false;
    @Parameter(names = "--ignore-input-index", description = "Require index of input", required = false)
    private boolean ignoreInputIndex = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        AddZgmId app = new AddZgmId();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private void go() {
        // TODO Auto-generated method stub
        String filename = inputVcfFile; // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        File file = new File(filename);

        VCFFileReader vcfFileReader = new VCFFileReader(file, false);
        VCFHeader header = vcfFileReader.getFileHeader();
        CloseableIterator<VariantContext> it = vcfFileReader.iterator();

        EnumSet<Options> options = ignoreInputIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet
                .noneOf(Options.class);
        if (ignoreMissingHeader)
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

        final VariantContextWriter vcfWriter = VariantContextWriterFactory.sortOnTheFly(
                VariantContextWriterFactory.create(new File(outputVcfFile), header.getSequenceDictionary(), options),
                10000);
        header.addMetaDataLine(new VCFInfoHeaderLine(ZGM_HASH_INFO_NAME, 1, VCFHeaderLineType.String,
                "unique identifier of this variant"));
        vcfWriter.writeHeader(header);

        int counter = 0;
        while (it.hasNext()) {
            VariantContext vc = it.next();

            vc.getCommonInfo().putAttribute(
                    ZGM_HASH_INFO_NAME,
                    VariantHashCounter.hash(vc.getChr(), vc.getStart(), vc.getEnd(), vc.getReference().getBaseString(),
                            vc.getAlternateAllele(0).getBaseString()));

            vcfWriter.add(vc);

            if (++counter % 10000 == 0) {
                log.debug("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfFileReader);
        vcfWriter.close();

    }
}
