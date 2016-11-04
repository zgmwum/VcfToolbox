package com.cloudinside.bio.VcfToolbox.annotators;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.cloudinside.bio.VcfToolbox.IAnnotator;
import com.cloudinside.bio.VcfToolbox.utils.ChromosomePosition;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

public class BedNonOverlappingAnnotator implements IAnnotator {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(BedNonOverlappingAnnotator.class);

    private final String columnName;
    private final String filename;

    private String infoPrefix;

    private final RangeMap<ChromosomePosition, String> map = TreeRangeMap.create();

    private final String description;

    /**
     * @param filename
     *            is a tabix preapred file
     */
    public BedNonOverlappingAnnotator(String filename, String columnName, String description) {

        this.description = description;
        this.columnName = columnName;
        this.filename = filename;

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                } else {
                    String[] splitted = line.split("\t");
                    String chr = splitted[0].intern();
                    String value = "";
                    Integer start = Integer.valueOf(splitted[1]);
                    Integer end = Integer.valueOf(splitted[2]);
                    if (splitted.length >= 4)
                        value = splitted[3].intern();

                    map.put(Range.openClosed(new ChromosomePosition(chr, start), new ChromosomePosition(chr, end)),
                            value);
                }
            }
            br.close();

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

    public Map<String, Object> getFeatures(String chr, int pos, String ref, String obs) {
        Map<String, Object> output = new HashMap<String, Object>();

        String value = map.get(new ChromosomePosition(chr, pos));
        if (value == null) {

        } else {
            output.put(columnName, value);
        }
        return output;
    }

    public Multimap<Integer, Map<String, Object>> getFeatures(String chr, int start, int stop) {

        Multimap<Integer, Map<String, Object>> multimap = HashMultimap.create();

        for (int i = start; i <= stop; i++) {
            String value = map.get(new ChromosomePosition(chr, i));
            if (value != null) {
                Map<String, Object> output = new HashMap<String, Object>();
                output.put(columnName, value);
                multimap.put(i, output);
            }
        }

        return multimap;
    }

    @Override
    public Collection<VCFInfoHeaderLine> getInfoHeaderLines() {

        Set<VCFInfoHeaderLine> lines = new HashSet<>();
        VCFInfoHeaderLine newHeaderLine = new VCFInfoHeaderLine(columnName, 1, VCFHeaderLineType.Character, description);
        lines.add(newHeaderLine);
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
                    log.debug("Duplicated key " + newKey);
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
                    log.debug("Duplicated key " + newKey);
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

    @Override
    public void close() throws IOException {
        // TODO Auto-generated method stub

    }

    @Override
    public void setSelectAttributes(Collection<String> attributes) {
        // TODO Auto-generated method stub

    }

    private String infoSuffix;

    @Override
    public void setSuffix(String infoSuffix) {
        this.infoSuffix = infoSuffix;

    }

    public static void main(String[] args) {
        BedNonOverlappingAnnotator ba = new BedNonOverlappingAnnotator(
                "/wum/pio/annotatios/mito-annotations/mtLoc.bed", "column", "desc");
        VariantContext vs = new VariantContextBuilder().chr("chrM").start(112).stop(112).alleles("A", "C").make();
        VariantContext avs = ba.annotate(vs);
        System.out.println(avs);
    }
}
