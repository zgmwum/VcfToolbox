package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.SortingVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Read multisample vcf file
 * 
 * Writes file with no samples, INFO contains sample names and frequencies of
 * observed allele. Multiallele obs are splitted
 * 
 * @author pstawinski
 * 
 */
public class BuildFrequencySummary {
    private static final String REFERENCE_NUMBER_INFO = "RefNum";
    private static final String ALLELE_FREQUENCY_INFO = "AllFreq";
    private static final String ALL_ALLELES_COUNT_INFO = "TotNum";
    private static final String SAMPLES_COUNT_INFO = "SampNum";
    private static final String SAMPLES_CONTAINING_INFO = "SampCont";
    private static final String SAMPLES_CONTAINING_FREQUENCY_INFO = "SampFreq";
    private static final String SAMPLES_INFO = "Samples";
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(BuildFrequencySummary.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--prefix", description = "Prefix of info entries", required = false)
    private String infoPrefix = "ZGM_";
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header", required = false)
    private boolean ignoreMissingHeader = false;
    @Parameter(names = "--normalize", description = "Normalize alleles, usually not needed", required = false)
    private boolean normalize = false;
    // @Parameter(names = "--ignore-input-index", description = "Require index
    // of input", required = false)
    private boolean ignoreInputIndex = true;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        BuildFrequencySummary app = new BuildFrequencySummary();
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

        String filename = inputVcfFile; // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        File file = new File(filename);

        VCFFileReader vcfFileReader = new VCFFileReader(file, !ignoreInputIndex);
        VCFHeader header = vcfFileReader.getFileHeader();
        CloseableIterator<VariantContext> it = vcfFileReader.iterator();

        final EnumSet<Options> options = ignoreInputIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY)
                : EnumSet.noneOf(Options.class);

        options.add(Options.USE_ASYNC_IO);

        final VariantContextWriter vcfWriter = new SortingVariantContextWriter(
                new VariantContextWriterBuilder().setOptions(options).setOutputFile(new File(outputVcfFile))
                        .setReferenceDictionary(header.getSequenceDictionary()).build(),
                10000);

        Set<VCFHeaderLine> set = Collections.emptySet();
        List<String> list = Collections.emptyList();
        VCFHeader vcfHeader = new VCFHeader(set, list);

        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + ALL_ALLELES_COUNT_INFO, 1,
                VCFHeaderLineType.Integer, "Total number of observations at this position"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + ALLELE_FREQUENCY_INFO, 1, VCFHeaderLineType.Float,
                "Current allele frequency divided by number of observations"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + REFERENCE_NUMBER_INFO, 1,
                VCFHeaderLineType.Integer, "Number of reference equals observations at this point"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_INFO, VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String, "Sample names, that included this allele"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_COUNT_INFO, 1, VCFHeaderLineType.Integer,
                "All samples that were analyzed"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_CONTAINING_INFO, 1,
                VCFHeaderLineType.Integer, "Samples containing this allele"));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_CONTAINING_FREQUENCY_INFO, 1,
                VCFHeaderLineType.Float, "Fraction of samples containing this allele"));

        // vcfHeader.setSequenceDictionary(header.getSequenceDictionary());

        // vcfHeader.addMetaDataLine()

        vcfWriter.writeHeader(vcfHeader);

        int counter = 0;
        while (it.hasNext()) {
            try {
                VariantContext vc = it.next();

                try {

                    Multimap<Allele, Genotype> alleles = ArrayListMultimap.create();
                    int referenceCount = 0;

                    // System.out.println(vc);
                    String chr = vc.getChr();
                    Allele oryginalReference = vc.getReference();
                    Set<String> sampleNames = vc.getSampleNames();
                    int allSamplesCount = sampleNames.size();

                    GenotypesContext genotypesContext = vc.getGenotypes();
                    for (Genotype genotype : genotypesContext) {
                        for (Allele allele : genotype.getAlleles()) {
                            if (allele.isReference()) {
                                referenceCount++;
                            } else {
                                if (allele.isCalled()) {
                                    alleles.put(allele, genotype);
                                } else {
                                    // no called, probably not available here
                                }
                            }
                        }
                    }

                    for (Entry<Allele, Collection<Genotype>> alleleToGenotypes : alleles.asMap().entrySet()) {
                        Collection<Genotype> genotypes = alleleToGenotypes.getValue();
                        Allele allele = alleleToGenotypes.getKey();

                        VariantContextBuilder vcb = new VariantContextBuilder();
                        vcb.chr(chr);
                        vcb.attributes(new HashMap<String, Object>());
                        vcb.attribute(infoPrefix + REFERENCE_NUMBER_INFO, referenceCount);
                        vcb.attribute(infoPrefix + ALL_ALLELES_COUNT_INFO, referenceCount + genotypes.size());

                        Set<String> samplesList = new TreeSet<>();
                        for (Genotype g : genotypes) {
                            samplesList.add(g.getSampleName());
                        }
                        vcb.attribute(infoPrefix + SAMPLES_INFO, Joiner.on(',').join(samplesList));
                        vcb.attribute(infoPrefix + SAMPLES_COUNT_INFO, allSamplesCount);
                        vcb.attribute(infoPrefix + SAMPLES_CONTAINING_INFO, samplesList.size());
                        vcb.attribute(infoPrefix + SAMPLES_CONTAINING_FREQUENCY_INFO,
                                (double) samplesList.size() / allSamplesCount);

                        vcb.attribute(infoPrefix + ALLELE_FREQUENCY_INFO,
                                (double) genotypes.size() / (referenceCount + genotypes.size()));

                        int position = vc.getStart();
                        String refSequence = oryginalReference.getBaseString();
                        String alleleSequence = allele.getBaseString();
                        if (normalize) {

                            int prefixOffset = countPrefixOffset(refSequence, alleleSequence);
                            refSequence = refSequence.substring(prefixOffset);
                            alleleSequence = alleleSequence.substring(prefixOffset);

                            int suffixOffset = countSuffixOffset(refSequence, alleleSequence);
                            refSequence = refSequence.substring(0, refSequence.length() - suffixOffset);
                            alleleSequence = alleleSequence.substring(0, alleleSequence.length() - suffixOffset);

                            position += prefixOffset;
                        }

                        vcb.start(position).alleles(refSequence, alleleSequence);
                        vcb.noGenotypes().noID();
                        vcb.computeEndFromAlleles(vcb.getAlleles(), position);

                        VariantContext outputVc = vcb.make();

                        vcfWriter.add(outputVc);
                    }
                } catch (Exception e) {
                    System.err.println("Error for: " + vc);
                    e.printStackTrace();
                }
            } catch (Exception e) {
                System.err.println("Ommiting line due to ar error");
                e.printStackTrace();
            }

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfFileReader);
        vcfWriter.close();

    }

    private int countPrefixOffset(String oryginalReference, String allele) {
        if (oryginalReference.length() == 1 || allele.length() == 1) {
            return 0;
        }

        String ref = oryginalReference;
        String obs = allele;
        String commonPrefix = Strings.commonPrefix(ref, obs);

        int offset = commonPrefix.length() - 1; // we want to leave last common
                                                // base at the beginning
        if (offset < 0) {
            return 0;
        } else {
            return offset;
        }

    }

    private int countSuffixOffset(String oryginalReference, String allele) {
        if (oryginalReference.length() == 1 || allele.length() == 1) {
            return 0;
        }
        String ref = oryginalReference;
        String obs = allele;
        String commonSuffix = Strings.commonSuffix(ref, obs);

        int offset = commonSuffix.length();
        if (offset == allele.length() || offset == oryginalReference.length()) {
            return offset - 1;
        } else {
            return offset;
        }
    }
}
