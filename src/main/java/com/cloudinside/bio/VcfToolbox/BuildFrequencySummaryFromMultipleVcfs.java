package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultimap;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.SortingVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Read multiple single sample trimmed, left aligned vcfs from stdin as text.
 * Writes vcf with appropriate information
 * 
 * @author pstawinski
 */
public class BuildFrequencySummaryFromMultipleVcfs {
    private static final String REFERENCE_NUMBER_INFO = "RefNum";
    private static final String ALLELE_FREQUENCY_INFO = "AllFreq";
    private static final String ALL_ALLELES_COUNT_INFO = "TotNum";
    private static final String SAMPLES_COUNT_INFO = "SampNumber";
    private static final String SAMPLES_CONTAINING_INFO = "SampContaining";
    private static final String SAMPLES_CONTAINING_FREQUENCY_INFO = "SampFreq";
    private static final String SAMPLES_INFO = "Samples";
    private static final String SAMPLES_PHENOTYPES = "Phen";
    private static final String SAMPLES_PHENOTYPES_HISTOGRAM = "PhenCount";
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(BuildFrequencySummaryFromMultipleVcfs.class);

    // @Parameter(names = "--sequence-dictionary-source", description =
    // "Vcf input file (can be bgzipped); is used as a source of reference
    // dictionary data",
    // required = true)
    // private String inputVcfFile;
    @Parameter(names = "--prefix", description = "Prefix of info entries", required = false)
    private String infoPrefix = "";
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header", required = false)
    private boolean ignoreMissingHeader = true;
    @Parameter(names = "--ignore-input-index", description = "Require index of input", required = false)
    private boolean ignoreInputIndex = false;
    @Parameter(names = "--strip-common-name-prefix", description = "Strip common prefix from sample names, f.ex. Sample_", required = false)
    private String commonNamePrefixToStrip = "";
    @Parameter(names = "--count-zgm-phenotype-freq", description = "Make a histogram of zgm-phenotype in variants", required = false)
    private boolean countZgmPhenotype = false;
    @Parameter(names = "--mark-homo", description = "Mark homozygotes with a *", required = false)
    private boolean markHomo = true;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        BuildFrequencySummaryFromMultipleVcfs app = new BuildFrequencySummaryFromMultipleVcfs();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    static NumberFormat df = DecimalFormat.getInstance(Locale.ENGLISH);
    static {
        df.setMinimumFractionDigits(0);
        df.setMaximumFractionDigits(4);
    }

    private void go() {

        // TODO Auto-generated method stub
        // String filename = inputVcfFile; //
        // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        // File file = new File(filename);

        // VCFFileReader vcfFileReader = new VCFFileReader(file,
        // !ignoreInputIndex);
        // VCFHeader header = vcfFileReader.getFileHeader();

        final EnumSet<Options> options = ignoreInputIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY)
                : EnumSet.noneOf(Options.class);
        // final VariantContextWriter vcfWriter = VariantContextWriterFactory
        // .sortOnTheFly(VariantContextWriterFactory.create(new
        // File(outputVcfFile), null, options), 1000);

        final VariantContextWriter vcfWriter = new SortingVariantContextWriter(
                new VariantContextWriterBuilder().setOutputFile(new File(outputVcfFile)).setOptions(options).build(),
                10000, true);

        Set<VCFHeaderLine> set = Collections.emptySet();
        List<String> list = Collections.emptyList();
        VCFHeader vcfHeader = new VCFHeader(set, list);
        // vcfHeader.setSequenceDictionary(header.getSequenceDictionary());

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

        if (countZgmPhenotype) {
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_PHENOTYPES, 1,
                    VCFHeaderLineType.String, "Phenotypes where this variant appears"));
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(infoPrefix + SAMPLES_PHENOTYPES_HISTOGRAM, 1,
                    VCFHeaderLineType.String, "Histogram of sample phenotypes"));
        }

        // vcfHeader.setSequenceDictionary(header.getSequenceDictionary());

        // vcfHeader.addMetaDataLine()

        vcfWriter.writeHeader(vcfHeader);

        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

        int counter = 0;
        int samplesCounter = 0;
        String line;
        String currentSampleName = null;
        String currentPhenotype = "NN";
        Map<String, String> sampleNameToPhenotype = new HashMap<>();
        Multimap<VariantBasicData, String> variantToSampleName = TreeMultimap.create();
        try {
            while ((line = br.readLine()) != null) {
                if (line.startsWith("##")) {
                    // omit, continues
                    if (countZgmPhenotype && line.startsWith("##zgm_phenotype=")) {
                        currentPhenotype = line.substring("##zgm_phenotype=".length());
                    }
                } else if (line.startsWith("#CHROM")) {
                    String[] splitted = line.split("\\s");
                    String sampleName = splitted[splitted.length - 1];
                    if ("FORMAT".endsWith(sampleName)) {
                        log.error("Unable to determine sample name from line: " + line);
                        return;
                    }

                    String sampleNameStripped = sampleName;
                    if (StringUtils.isNotBlank(commonNamePrefixToStrip)) {
                        if (sampleName.startsWith(commonNamePrefixToStrip)) {
                            sampleNameStripped = sampleName.substring(commonNamePrefixToStrip.length());
                        }
                    }

                    currentSampleName = sampleNameStripped;

                    if (countZgmPhenotype)
                        sampleNameToPhenotype.put(currentSampleName, currentPhenotype);

                    samplesCounter++;

                    log.debug("Analyzing " + sampleName);

                } else {
                    String[] splitted = line.split("\\s");
                    String chr = splitted[0];
                    int pos = Integer.valueOf(splitted[1]);
                    // DBSNP [2]
                    String ref = splitted[3];
                    String obs = splitted[4];
                    // QUAL [5]
                    // FILTER [6]
                    // INFO [7]
                    String info = splitted[7];

                    Boolean homo = null;
                    if (markHomo) {
                        Map<String, String> infoMap = infoToStringMap(info);

                        if (infoMap.containsKey("AF")) {
                            Double alleleFrequency;
                            alleleFrequency = Double.valueOf(infoMap.get("AF"));
                            if (alleleFrequency.equals(1.)) {
                                homo = Boolean.TRUE;
                            } else if (alleleFrequency.equals(0.5)) {
                                homo = Boolean.FALSE;
                            }
                        }
                    }

                    VariantBasicData vbd = new VariantBasicData(chr, pos, ref, obs);
                    // check if is homo or heterozyg

                    variantToSampleName.put(vbd, currentSampleName + ((homo != null && homo) ? "*" : ""));
                }
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        log.debug("All samples loaded, saving vcf file");
        for (Entry<VariantBasicData, Collection<String>> entry : variantToSampleName.asMap().entrySet()) {
            VariantBasicData vbd = entry.getKey();
            List<String> sampleNames = new ArrayList<String>(entry.getValue());
            Collections.sort(sampleNames);

            VariantContextBuilder vcb = new VariantContextBuilder();
            vcb.chr(vbd.getChr()).start(vbd.getPos()).alleles(vbd.getRef(), vbd.getObs());
            vcb.attributes(new HashMap<String, Object>());

            vcb.attribute(infoPrefix + SAMPLES_INFO, Joiner.on(',').join(sampleNames));
            vcb.attribute(infoPrefix + SAMPLES_COUNT_INFO, samplesCounter);
            vcb.attribute(infoPrefix + SAMPLES_CONTAINING_INFO, sampleNames.size());
            vcb.attribute(infoPrefix + SAMPLES_CONTAINING_FREQUENCY_INFO,
                    df.format((double) sampleNames.size() / samplesCounter));

            if (countZgmPhenotype) {

                Multiset<String> phenotypes = HashMultiset.create();
                for (String sampleName : sampleNames) {
                    phenotypes.add(sampleNameToPhenotype.get(sampleName.replace("*", "")));
                }

                vcb.attribute(infoPrefix + SAMPLES_PHENOTYPES, Joiner.on(',').join(phenotypes.elementSet()));
                vcb.attribute(infoPrefix + SAMPLES_PHENOTYPES_HISTOGRAM, phenotypes.toString().replace(' ', '_'));

            }

            vcb.noGenotypes().noID();
            vcb.computeEndFromAlleles(vcb.getAlleles(), vbd.getPos());

            VariantContext outputVc = vcb.make();
            vcfWriter.add(outputVc);
        }

        IOUtils.closeQuietly(br);
        vcfWriter.close();

    }

    private Map<String, String> infoToStringMap(String info) {
        Map<String, String> map = new HashMap<>();
        for (String part : info.split(";")) {
            String[] partSplitted = part.split("=");
            if (partSplitted.length == 1) {
                map.put(partSplitted[0], "");
            } else {
                map.put(partSplitted[0], partSplitted[1]);
            }
        }
        return map;
    }
}
