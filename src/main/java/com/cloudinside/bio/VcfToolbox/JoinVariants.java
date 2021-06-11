package com.cloudinside.bio.VcfToolbox;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Vcf files must be in single obs form (no observations X,Y) sorted by
 * chromosome(alphanumeric)/position/ref/obs
 * 
 * <pre>
 * sort -k1,1V -k2,2n -k4,4 -k5,5
 * </pre>
 * 
 * @author pstawinski
 * 
 */
public class JoinVariants {
    private static final String ZGM_HASH_INFO_NAME = "ZgmHash";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(JoinVariants.class);

    private static final String ZGM_VCF_SUPPORTING_ANALYZED = "ZGM_vcf_sc";
    private static final String ZGM_VCF_SUPPORTING_FREQ = "ZGM_vcf_sf";
    private static final String ZGM_VCF_SUPPORTING_NUM = "ZGM_vcf_sn";
    private static final String ZGM_VCF_SUPPORTING_SOURCES = "ZGM_vcf_sources";

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private List<String> inputVcfFiles;
    @Parameter(names = "--prefixes", description = "Prefixes of vcfs, must be equal in size as --input", required = true)
    private List<String> inputPrefixes;
    @Parameter(names = "--ref", description = "Reference, used only for sequence dictionary", required = true)
    private String referenceFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header", required = false)
    private boolean ignoreMissingHeader = false;
    @Parameter(names = "--create-output-index", description = "Create index of output vcf file", required = false)
    private boolean createOutputIndex = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        JoinVariants app = new JoinVariants();
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
        if (inputPrefixes.size() != inputVcfFiles.size()) {
            System.err.println("inputPrefixes size must be equal inputVcfFiles size");
            throw new NullPointerException(); // :D
        }

        EnumSet<Options> options = createOutputIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY)
                : EnumSet.noneOf(Options.class);
        if (ignoreMissingHeader)
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriter vcfWriter = null;

        VCFHeader oputputHeader = null;

        List<VcfIteratorWrapper> readerIterators = new ArrayList<>();
        {
            int index = 0;
            for (String filename : inputVcfFiles) {
                File file = new File(filename);

                VCFFileReader vcfFileReader = new VCFFileReader(file, false);
                VCFHeader header = vcfFileReader.getFileHeader();
                CloseableIterator<VariantContext> it = vcfFileReader.iterator();
                readerIterators.add(new VcfIteratorWrapper(it));

                if (vcfWriter == null) {
                    vcfWriter = VariantContextWriterFactory.create(new File(outputVcfFile),
                            header.getSequenceDictionary(), options);

                }

                if (index == 0) {

                    Set<VCFHeaderLine> headerLines = new HashSet<>();
                    headerLines.addAll(header.getContigLines());
                    headerLines.addAll(header.getFilterLines());
                    headerLines.addAll(header.getFormatHeaderLines());
                    oputputHeader = new VCFHeader(headerLines,
                            Collections.singletonList(header.getSampleNamesInOrder().get(0)));

                    oputputHeader.addMetaDataLine(new VCFInfoHeaderLine(ZGM_VCF_SUPPORTING_NUM, 1,
                            VCFHeaderLineType.Integer, "number of files supporting this change"));
                    oputputHeader.addMetaDataLine(new VCFInfoHeaderLine(ZGM_VCF_SUPPORTING_FREQ, 1,
                            VCFHeaderLineType.Float, "frequency of files supporting this"));
                    oputputHeader.addMetaDataLine(new VCFInfoHeaderLine(ZGM_VCF_SUPPORTING_ANALYZED, 1,
                            VCFHeaderLineType.Integer, "number of vcf files analyzed"));
                    oputputHeader.addMetaDataLine(
                            new VCFInfoHeaderLine(ZGM_VCF_SUPPORTING_SOURCES, VCFHeaderLineCount.UNBOUNDED,
                                    VCFHeaderLineType.String, "vcfs supporting this observation"));
                }

                for (VCFInfoHeaderLine infoHeaderLine : header.getInfoHeaderLines()) {
                    VCFInfoHeaderLine newHeaderLine;
                    if (infoHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER)
                        newHeaderLine = new VCFInfoHeaderLine(inputPrefixes.get(index) + infoHeaderLine.getID(),
                                infoHeaderLine.getCount(), infoHeaderLine.getType(), infoHeaderLine.getDescription());
                    else
                        newHeaderLine = new VCFInfoHeaderLine(inputPrefixes.get(index) + infoHeaderLine.getID(),
                                infoHeaderLine.getCountType(), infoHeaderLine.getType(),
                                infoHeaderLine.getDescription());

                    oputputHeader.addMetaDataLine(newHeaderLine);
                }

                // merge Filter and Format headers
                for (VCFFilterHeaderLine filterHeaderLine : header.getFilterLines()) {
                    if (!oputputHeader.hasFilterLine(filterHeaderLine.getID())) {
                        oputputHeader.addMetaDataLine(filterHeaderLine);
                    }
                }
                for (VCFFormatHeaderLine formatHeaderLine : header.getFormatHeaderLines()) {
                    if (!oputputHeader.hasFormatLine(formatHeaderLine.getID())) {
                        oputputHeader.addMetaDataLine(formatHeaderLine);
                    }
                }

                // if (index == 0) {
                // oputputHeader.setSequenceDictionary(header.getSequenceDictionary());
                // }
                vcfFileReader.close();
                index++;
            }
        }

        vcfWriter.writeHeader(oputputHeader);

        ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
                .getReferenceSequenceFile(new File(referenceFile));

        VariantContextComparator variantContextComparator = new VariantContextComparator(
                referenceSequenceFile.getSequenceDictionary());
        try {
            referenceSequenceFile.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        int counter = 0;
        boolean done = false;
        while (!done) {
            // find smallest
            VariantContext smallest = null;
            for (VcfIteratorWrapper it : readerIterators) {
                if (smallest == null) {
                    smallest = it.peek();
                } else {
                    if (it.peek() != null && variantContextComparator.compare(smallest, it.peek()) > 0) {
                        smallest = it.peek();
                    }
                }
            }

            // smallest contains smallest or null
            if (smallest == null) {
                done = true;
            } else {
                VariantContext vc = smallest;
                VariantContextBuilder vcb = new VariantContextBuilder(vc);

                // vcb.filter(".");
                // vcb.noID();
                // vcb.noGenotypes();

                vcb.attributes(new HashMap<String, Object>());

                int count = 0;
                int index = 0;
                List<String> prefixes = new ArrayList<>();
                for (VcfIteratorWrapper it : readerIterators) {
                    if (it.peek() != null && variantContextComparator.compare(smallest, it.peek()) == 0) {
                        VariantContext itvc = it.pop();
                        for (Entry<String, Object> attribute : itvc.getAttributes().entrySet()) {
                            vcb.attribute(inputPrefixes.get(index) + attribute.getKey(), attribute.getValue());
                        }
                        count++;
                        prefixes.add(inputPrefixes.get(index));
                    }

                    index++;

                }
                Genotype genotype = vc.getGenotypes().get(0);
                vcb.genotypes(genotype);
                vcb.attribute(ZGM_VCF_SUPPORTING_ANALYZED, inputPrefixes.size());
                vcb.attribute(ZGM_VCF_SUPPORTING_NUM, count);
                vcb.attribute(ZGM_VCF_SUPPORTING_FREQ, (double) count / (double) inputPrefixes.size());
                vcb.attribute(ZGM_VCF_SUPPORTING_SOURCES, StringUtils.join(prefixes, ','));

                VariantContext outputVc = vcb.make();

                vcfWriter.add(outputVc);

            }

            // VariantContext vc = it.next();
            // Allele oryginalReference = vc.getReference();
            //
            // for (Allele obs : vc.getAlternateAlleles()) {
            // VariantContextBuilder vcb = new VariantContextBuilder(vc);
            //
            // vc.getCommonInfo().putAttribute(
            // ZGM_HASH_INFO_NAME,
            // VariantHashCounter.hash(vc.getChr(), vc.getStart(), vc.getEnd(),
            // vc.getReference()
            // .getBaseString(), vc.getAlternateAllele(0).getBaseString()));
            //
            // int position = vc.getStart();
            // String refSequence = oryginalReference.getBaseString();
            // String alleleSequence = obs.getBaseString();
            // int prefixOffset = countPrefixOffset(refSequence,
            // alleleSequence);
            // refSequence = refSequence.substring(prefixOffset);
            // alleleSequence = alleleSequence.substring(prefixOffset);
            // int suffixOffset = countSuffixOffset(refSequence,
            // alleleSequence);
            // refSequence = refSequence.substring(0, refSequence.length() -
            // suffixOffset);
            // alleleSequence = alleleSequence.substring(0,
            // alleleSequence.length() - suffixOffset);
            //
            // position += prefixOffset;
            //
            // vcb.start(position).alleles(refSequence, alleleSequence);
            // vcb.computeEndFromAlleles(vcb.getAlleles(), position);
            //
            // VariantContext outputVc = vcb.make();
            // vcfWriter.add(outputVc);
            // }

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        for (VcfIteratorWrapper it : readerIterators) {
            IOUtils.closeQuietly(it);
        }
        vcfWriter.close();

    }

    private static class VcfIteratorWrapper implements Closeable {
        private final CloseableIterator<VariantContext> iterator;

        private VariantContext topElement;

        public VcfIteratorWrapper(CloseableIterator<VariantContext> iterator) {
            this.iterator = iterator;
        }

        public VariantContext peek() {
            if (topElement != null)
                return topElement;
            else
                loadElementToTopElement();

            return topElement;
        }

        public VariantContext pop() {
            if (topElement != null) {
                VariantContext vc = topElement;
                topElement = null;
                return vc;
            } else {
                loadElementToTopElement();
                if (topElement != null) {
                    VariantContext vc = topElement;
                    topElement = null;
                    return vc;
                } else {
                    return null;
                }
            }
        }

        private void loadElementToTopElement() {
            if (iterator.hasNext())
                topElement = iterator.next();
            else
                topElement = null;
        }

        @Override
        public void close() throws IOException {
            iterator.close();
        }
    }

}
