package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class CopyAttributesFromFirstVariantToSameVariants {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(CopyAttributesFromFirstVariantToSameVariants.class);

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

        CopyAttributesFromFirstVariantToSameVariants app = new CopyAttributesFromFirstVariantToSameVariants();
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

        vcfWriter.writeHeader(header);

        String prevRef = null, prevObs = null, prevChr = null;
        int prevPos = -1;
        VariantContext prevVariantContext = null;

        int counter = 0;
        while (it.hasNext()) {
            VariantContext vc = it.next();
            String curRef = vc.getReference().getBaseString();
            String curObs = Joiner.on(',').join(vc.getAlternateAlleles());
            String curChr = vc.getChr();
            int curPos = vc.getStart();

            if (prevVariantContext != null) {
                if (prevPos == curPos && curRef.equalsIgnoreCase(prevRef) && curObs.equalsIgnoreCase(prevObs)
                        && curChr.equalsIgnoreCase(prevChr)) {
                    VariantContextBuilder vcb = new VariantContextBuilder(vc);
                    Map<String, Object> currentAttributes = vc.getCommonInfo().getAttributes();
                    for (Entry<String, Object> prevAttribute : prevVariantContext.getAttributes().entrySet()) {
                        if (currentAttributes.containsKey(prevAttribute.getKey())) {
                            // omit
                        } else {
                            vc.getCommonInfo().putAttribute(prevAttribute.getKey(), prevAttribute.getValue());
                        }
                    }
                    if (vc.emptyID()) {
                        vcb.id(prevVariantContext.getID());
                    }
                    VariantContext outputVc = vcb.make();
                    vcfWriter.add(outputVc);
                    vc = outputVc;

                } else {
                    vcfWriter.add(vc);

                }
            }
            prevVariantContext = vc;

            prevRef = curRef;
            prevObs = curObs;
            prevPos = curPos;
            prevChr = curChr;

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfFileReader);
        vcfWriter.close();

    }
}
