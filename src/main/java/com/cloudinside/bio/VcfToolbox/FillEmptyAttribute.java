package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;

import net.sf.samtools.util.CloseableIterator;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Stopwatch;

/**
 * Finds variants in genes that are possibly damaging and appears twice (i.e.
 * possibly two gene copies are destroyed).
 * 
 * @author pstawinski
 * 
 */
public class FillEmptyAttribute {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(FillEmptyAttribute.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix). Cann vcf file with snpEff annotations", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header lines", required = false)
    private boolean ignoreMissingHeader = false;

    @Parameter(names = "--attribute", description = "Attribute name", required = true)
    private List<String> attributes;
    @Parameter(names = "--value", description = "Attribute value", required = true)
    private List<String> values;

    // @Parameter(names = "--ignore-input-index", description =
    // "Require index of input", required = false)
    private boolean ignoreInputIndex = true;

    private CloseableIterator<VariantContext> it;

    private Stopwatch stopwatch;

    private int counter;

    private VariantContextWriter vcfWriter;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        FillEmptyAttribute app = new FillEmptyAttribute();
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

    private void go() {
        // TODO Auto-generated method stub
        String filename = inputVcfFile; // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        File file = new File(filename);

        VCFFileReader vcfFileReader = new VCFFileReader(file, false);
        VCFHeader header = vcfFileReader.getFileHeader();
        it = vcfFileReader.iterator();

        EnumSet<Options> options = EnumSet.noneOf(Options.class);
        if (ignoreMissingHeader)
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

        vcfWriter = // VariantContextWriterFactory.sortOnTheFly(
        VariantContextWriterFactory.create(new File(outputVcfFile), header.getSequenceDictionary(), options);
        // 10000);

        // build header

        vcfWriter.writeHeader(header);

        counter = 0;
        stopwatch = Stopwatch.createStarted();

        it = vcfFileReader.iterator();

        {
            VariantContext vc;
            while (it.hasNext()) {
                vc = it.next();
                VariantContextBuilder vcb = new VariantContextBuilder(vc);

                int index = 0;
                for (String attribute : attributes) {
                    Object attributeValue = vc.getAttribute(attribute);
                    if (attributeValue == null) {
                        vcb.attribute(attribute, values.get(index));
                    }
                    index++;
                }

                vcfWriter.add(vcb.make());

            }
        }

        vcfWriter.close();

    }
}
