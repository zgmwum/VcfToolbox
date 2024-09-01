package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.Locale;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Read single or multisample vcf file. Remove these rows, where number of read
 * for ANY of the sample is too low.
 * 
 * @author pstawinski
 * 
 */
public class RemoveLowCoverageRows {

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--min-dp", description = "Minimal DP", required = false)
    private int minDp = 10;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        RemoveLowCoverageRows app = new RemoveLowCoverageRows();
        JCommander jc = null;
        try {
            jc = new JCommander(app);
            jc.parse(args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private void go() {
        File vcfFile = new File(inputVcfFile);
        File outputFile = new File(outputVcfFile);

        final EnumSet<Options> options = EnumSet.noneOf(Options.class);

        options.add(Options.USE_ASYNC_IO);

        try (VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
                VariantContextWriter vcfWriter = new VariantContextWriterBuilder().setOptions(options)
                        .setOutputFile(outputFile)

                        .build()) {

            vcfWriter.writeHeader(vcfReader.getFileHeader());

            int counter = 0;

            for (VariantContext vc : vcfReader) {
                boolean keepVariant = false;

                // Check DP for each sample
                for (Genotype genotype : vc.getGenotypes()) {
                    int dp = 0;
                    if (genotype.hasDP()) {
                        dp = genotype.getDP();
                    }

                    // If any sample has DP >= dpThreshold, keep the variant
                    if (dp >= minDp) {
                        keepVariant = true;
                        break;
                    }
                }

                if (keepVariant) {
                    vcfWriter.add(vc);
                }

                if (++counter % 10000 == 0) {
                    System.err.println("Processed " + counter);
                }
            }

        } catch (Exception e) {
            System.err.println("Error processing VCF file");
            e.printStackTrace();
        }

    }

}
