package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Read multisample vcf file, writes multiple vcfs, one per sample
 * 
 * @author pstawinski
 * 
 */
public class SplitMultisampleVcf {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(SplitMultisampleVcf.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output-prefix", description = "Prefix to be added to vcf outputs, for example /tmp/", required = false)
    private String outputPrefix = "";

    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header", required = false)
    private boolean ignoreMissingHeader = false;
    @Parameter(names = "--ignore-input-index", description = "Require index of input", required = false)
    private boolean ignoreInputIndex = true;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        SplitMultisampleVcf app = new SplitMultisampleVcf();
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

        List<String> sampleNames = vcfFileReader.getFileHeader().getSampleNamesInOrder();
        Map<String, VariantContextWriter> writers = new HashMap<>();
        for (String sampleName : sampleNames) {
            final VariantContextWriter vcfWriter = VariantContextWriterFactory.sortOnTheFly(VariantContextWriterFactory
                    .create(new File(outputPrefix + sampleName + ".vcf"), header.getSequenceDictionary(), options),
                    10000);
            vcfWriter.writeHeader(header);
            writers.put(sampleName, vcfWriter);

        }

        int counter = 0;
        while (it.hasNext()) {
            VariantContext vc = it.next();
            for (String sampleName : vc.getSampleNames()) {

                Genotype genotype = vc.getGenotype(sampleName);
                if (genotype.isCalled()) {
                    VariantContextBuilder vcb = new VariantContextBuilder(vc);
                    vcb.genotypes(Collections.singletonList(genotype));
                    writers.get(sampleName).add(vcb.make());
                }

            }

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfFileReader);
        for (VariantContextWriter vcw : writers.values()) {
            vcw.close();
        }

    }

}
