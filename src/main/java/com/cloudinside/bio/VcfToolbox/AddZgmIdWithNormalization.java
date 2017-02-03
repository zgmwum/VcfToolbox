package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.Locale;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VariantHashCounter;
import com.google.common.base.Strings;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class AddZgmIdWithNormalization {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(AddZgmIdWithNormalization.class);

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

        AddZgmIdWithNormalization app = new AddZgmIdWithNormalization();
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

        EnumSet<Options> options = ignoreInputIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY)
                : EnumSet.noneOf(Options.class);
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
            Allele oryginalReference = vc.getReference();

            for (Allele obs : vc.getAlternateAlleles()) {
                VariantContextBuilder vcb = new VariantContextBuilder(vc);

                vc.getCommonInfo().putAttribute(ZGM_HASH_INFO_NAME, VariantHashCounter.hash(vc.getChr(), vc.getStart(),
                        vc.getEnd(), vc.getReference().getBaseString(), vc.getAlternateAllele(0).getBaseString()));

                int position = vc.getStart();
                String refSequence = oryginalReference.getBaseString();
                String alleleSequence = obs.getBaseString();
                int prefixOffset = countPrefixOffset(refSequence, alleleSequence);
                refSequence = refSequence.substring(prefixOffset);
                alleleSequence = alleleSequence.substring(prefixOffset);
                int suffixOffset = countSuffixOffset(refSequence, alleleSequence);
                refSequence = refSequence.substring(0, refSequence.length() - suffixOffset);
                alleleSequence = alleleSequence.substring(0, alleleSequence.length() - suffixOffset);

                position += prefixOffset;

                vcb.start(position).alleles(refSequence, alleleSequence);
                vcb.computeEndFromAlleles(vcb.getAlleles(), position);

                VariantContext outputVc = vcb.make();
                vcfWriter.add(outputVc);
            }

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfFileReader);
        vcfWriter.close();

    }

    private int countPrefixOffset(String oryginalReference, String allele) {
        String ref = oryginalReference;
        String obs = allele;
        String commonPrefix = Strings.commonPrefix(ref, obs);

        int offset = commonPrefix.length() - 1;
        if (offset < 0) {
            return 0;
        } else {
            return offset;
        }

    }

    private int countSuffixOffset(String oryginalReference, String allele) {
        String ref = oryginalReference;
        String obs = allele;
        String commonSuffix = Strings.commonSuffix(ref, obs);

        int offset = commonSuffix.length() - 1;
        if (offset < 0) {
            return 0;
        } else {
            return offset;
        }
    }
}
