package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;

import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

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
public class RecessiveMarker {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(RecessiveMarker.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix). Cann vcf file with snpEff annotations", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header lines", required = false)
    private boolean ignoreMissingHeader = false;
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

        RecessiveMarker app = new RecessiveMarker();
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

        Set<String> transcriptWithImportantVariantSeen = new HashSet<>();
        Set<String> transcriptWithImportantVariantSeenTwice = new HashSet<>();

        vcfWriter = // VariantContextWriterFactory.sortOnTheFly(
        VariantContextWriterFactory.create(new File(outputVcfFile), header.getSequenceDictionary(), options);
        // 10000);

        // build header

        header.addMetaDataLine(new VCFInfoHeaderLine("ZGM_PROB_REC", 1, VCFHeaderLineType.Flag,
                "Variant in gene with at least two important variants"));

        vcfWriter.writeHeader(header);

        counter = 0;
        stopwatch = Stopwatch.createStarted();

        {
            VariantContext vc;
            while (it.hasNext()) {
                vc = it.next();
                Object effObj = vc.getAttribute("EFF");
                String eff = (String) effObj;

                if (StringUtils.isBlank(eff))
                    continue;
                eff = eff.substring(4);
                // String topEff = eff.substring(0, eff.indexOf('('));
                String[] effAttr = eff.substring(eff.indexOf('(') + 1, eff.length() - 2).split("\\|");
                // String impact = effAttr[0];
                if (effAttr.length > TRANSCRIPT_NAME_POSITION) {
                    String geneName = effAttr[TRANSCRIPT_NAME_POSITION];
                    // if ("MODERATE".equals(impact) || "HIGH".equals(impact) ||
                    // "MODIFIER".equals(impact)) {
                    if (transcriptWithImportantVariantSeen.add(geneName)) {
                        transcriptWithImportantVariantSeenTwice.add(geneName);
                    }
                    // }
                }

            }
        }

        it.close();

        // SECOND PASS

        // IOUtils.closeQuietly(vcfFileReader);
        // vcfFileReader = new VCFFileReader(file, false);
        it = vcfFileReader.iterator();

        {
            VariantContext vc;
            while (it.hasNext()) {
                vc = it.next();
                String eff = (String) vc.getAttribute("EFF");
                if (StringUtils.isBlank(eff))
                    continue;
                eff = eff.substring(4);
                // String topEff = eff.substring(0, eff.indexOf('('));
                String[] effAttr = eff.substring(eff.indexOf('(') + 1, eff.length() - 2).split("\\|");
                // String impact = effAttr[0];
                boolean add = false;

                if (effAttr.length > TRANSCRIPT_NAME_POSITION) {
                    String geneName = effAttr[TRANSCRIPT_NAME_POSITION];
                    // if ("MODERATE".equals(impact) || "HIGH".equals(impact) ||
                    // "MODIFIER".equals(impact)) {
                    if (transcriptWithImportantVariantSeenTwice.contains(geneName)) {
                        add = true;
                    } else if (vc.getAttributeAsBoolean("HOM", false)) {
                        add = true;
                    }
                    // }
                }
                VariantContextBuilder vcb = new VariantContextBuilder(vc);
                if (add) {
                    vcb.attribute("ZGM_PROB_REC", true);
                }

                vcfWriter.add(vcb.make());

            }
        }
        IOUtils.closeQuietly(vcfFileReader);

        vcfWriter.close();

    }

}
