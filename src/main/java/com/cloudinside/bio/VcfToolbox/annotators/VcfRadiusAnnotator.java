package com.cloudinside.bio.VcfToolbox.annotators;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.cloudinside.bio.VcfToolbox.IAnnotator;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.utils.RangeMultimap;
import com.google.common.base.Joiner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

public class VcfRadiusAnnotator implements IAnnotator {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(VcfRadiusAnnotator.class);

    // private final Collection<String> columnNames;
    private final String filename;

    private final Collection<String> columnNames = null;

    private String infoPrefix = "";

    private final Map<String, RangeMultimap<LinkedHashMap<String, Object>>> chrToAnnotationsMultimap = new HashMap<>();

    private int radius;

    private VCFHeader header;

    private boolean stripChrFromChromosomeName;

    /**
     * @param filename
     *            is a tabix preapred file
     */
    public VcfRadiusAnnotator(String filename, boolean stripChrFromChromosomeName, int radius) {

        // this.columnNames = columnNames;
        this.filename = filename;
        this.radius = radius;
        this.stripChrFromChromosomeName = stripChrFromChromosomeName;

        try {
            BufferedReader br;
            if (filename.endsWith(".gz"))
                br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
            else
                br = new BufferedReader(new FileReader(filename));
            VcfStreamReader vcfStreamReader = new VcfStreamReader(br);
            header = vcfStreamReader.getVcfHeader();

            VcfLine line;
            while ((line = vcfStreamReader.next()) != null) {
                String chr = line.getChr();
                if (stripChrFromChromosomeName) {
                    chr = StringUtils.removeStart(chr, "chr");
                } else if (!chr.startsWith("chr")) {
                    chr = "chr" + chr;

                }
                chr = chr.intern();

                String value = "";
                Integer start = line.getPosition() - radius;
                Integer end = line.getPosition() + radius;

                RangeMultimap<LinkedHashMap<String, Object>> annotationsMultimap = chrToAnnotationsMultimap.get(chr);
                if (annotationsMultimap == null) {
                    annotationsMultimap = new RangeMultimap<>();
                    chrToAnnotationsMultimap.put(chr, annotationsMultimap);
                }

                annotationsMultimap.put(Range.closed(start, end), line.getInfo());

            }
            br.close();
            vcfStreamReader.close();

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    @Override
    public Collection<Map<String, Object>> getPointData(String chr, int pos) {
        return getFeatures(chr, pos, pos).get(pos);
    }

    private final Joiner joiner = Joiner.on(",");

    public Map<String, Object> getFeatures(String chr, int pos, String ref, String obs) {
        Map<String, Object> output = new HashMap<String, Object>();

        if (stripChrFromChromosomeName) {
            chr = StringUtils.removeStart(chr, "chr");
        }

        Set<LinkedHashMap<String, Object>> value = null;
        RangeMultimap<LinkedHashMap<String, Object>> multimap = chrToAnnotationsMultimap.get(chr);
        if (multimap != null && multimap.get(pos) != null && !multimap.get(pos).isEmpty())
            value = multimap.get(pos);
        if (value == null) {

        } else {
            Map<String, Collection<Object>> collapsed = collapse(value);

            for (Entry<String, Collection<Object>> entry : collapsed.entrySet()) {
                output.put(entry.getKey(), joiner.join(entry.getValue()));
            }
        }
        return output;
    }

    private Map<String, Collection<Object>> collapse(Set<LinkedHashMap<String, Object>> value) {
        Multimap<String, Object> multimap = ArrayListMultimap.create();
        for (LinkedHashMap<String, Object> part : value) {
            for (Entry<String, Object> entry : part.entrySet()) {
                multimap.put(entry.getKey(), entry.getValue());
            }
        }

        return multimap.asMap();
    }

    public Multimap<Integer, Map<String, Object>> getFeatures(String chr, int start, int stop) {

        throw new NullPointerException();
    }

    @Override
    public Collection<VCFInfoHeaderLine> getInfoHeaderLines() {

        VCFFileReader vcfFileReader = new VCFFileReader(new File(filename), false);
        VCFHeader headerLocal = vcfFileReader.getFileHeader();
        vcfFileReader.close();

        Set<VCFInfoHeaderLine> lines = new HashSet<>();
        for (VCFInfoHeaderLine infoHeaderLine : headerLocal.getInfoHeaderLines()) {
            VCFInfoHeaderLine newHeaderLine;
            if (infoHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER)
                newHeaderLine = new VCFInfoHeaderLine(infoPrefix + infoHeaderLine.getID() + infoSuffix,
                        infoHeaderLine.getCount(), infoHeaderLine.getType(), infoHeaderLine.getDescription());
            else
                newHeaderLine = new VCFInfoHeaderLine(infoPrefix + infoHeaderLine.getID() + infoSuffix,
                        infoHeaderLine.getCountType(), infoHeaderLine.getType(), infoHeaderLine.getDescription());

            if (columnNames == null || columnNames.isEmpty() || columnNames.contains(infoHeaderLine.getID())) {
                lines.add(newHeaderLine);
            }
        }

        return lines;
    }

    @Override
    public VariantContext annotate(VariantContext vc) {
        VariantContextBuilder vcb = new VariantContextBuilder(vc);
        Map<String, Object> attributes = vc.getAttributes();

        for (Allele allele : vc.getAlternateAlleles()) {
            Map<String, Object> features = getFeatures(vc.getChr(), vc.getStart(), vc.getReference().getBaseString(),
                    allele.getBaseString());

            for (Entry<String, Object> entry : features.entrySet()) {
                String key = entry.getKey();

                String newKey = infoPrefix + key + infoSuffix;
                if (!attributes.containsKey(newKey)) {
                    vcb.attribute(newKey, entry.getValue());
                } else {
                    vcb.attribute(newKey, entry.getValue());
                    vcb.attribute(newKey, attributes.get(newKey) + ";" + entry.getValue());
                    // log.debug("Duplicated key " + newKey);
                }
            }
        }

        return vcb.make();
    }

    @Override
    public VcfLine annotate(VcfLine vc) {

        for (String allele : vc.getAlt()) {
            Map<String, Object> features = getFeatures(vc.getChr(), vc.getPosition(), vc.getRef(), allele);

            for (Entry<String, Object> entry : features.entrySet()) {
                String key = entry.getKey();

                String newKey = infoPrefix + key + infoSuffix;
                if (!vc.getInfo().containsKey(newKey)) {
                    vc.getInfo().put(newKey, entry.getValue());
                } else {
                    vc.getInfo().put(newKey, vc.getInfo().get(newKey) + ";" + entry.getValue());
                    // log.debug("Duplicated key " + newKey);
                }
            }
        }

        return vc;
    }

    @Override
    public Collection<VariantContext> annotateMultiplePossible(Collection<VariantContext> vcs) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void setPrefix(String infoPrefix) {
        this.infoPrefix = infoPrefix;

    }

    private String infoSuffix = "";

    @Override
    public void setSuffix(String infoSuffix) {
        this.infoSuffix = infoSuffix;

    }

    @Override
    public void close() throws IOException {
        // TODO Auto-generated method stub

    }

    @Override
    public void setSelectAttributes(Collection<String> attributes) {
        // TODO Auto-generated method stub

    }

    public static void main(String[] args) throws IOException {
        VcfRadiusAnnotator ba = new VcfRadiusAnnotator(
                "/wum/pio/annotatios/HGMD/2014.1/Genome_Trax_hg19_vcf/hgmd_hg19.sort.vcf.gz", false, 10);
        VariantContext vs = new VariantContextBuilder().chr("chr1").start(874814).stop(874814).alleles("A", "C").make();
        VariantContext avs = ba.annotate(vs);
        ba.close();
        System.out.println(avs);
    }
}
