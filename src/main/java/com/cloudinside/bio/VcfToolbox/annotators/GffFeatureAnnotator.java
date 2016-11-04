package com.cloudinside.bio.VcfToolbox.annotators;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.broad.tribble.readers.TabixReader.Iterator;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.cloudinside.bio.VcfToolbox.IAnnotator;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.google.common.base.Splitter;
import com.google.common.base.Splitter.MapSplitter;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

public class GffFeatureAnnotator extends TabixAnnotator implements IAnnotator {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(GffFeatureAnnotator.class);

    private final Collection<String> columnNames;
    private final String filename;
    private final MapSplitter mapSplitter = Splitter.onPattern("; *").withKeyValueSeparator('=');

    private String infoPrefix = "";

    /**
     * @param filename
     *            is a tabix preapred file
     */
    public GffFeatureAnnotator(String filename, Collection<String> columnNames) {
        super(filename);
        this.columnNames = columnNames;
        this.filename = filename;
    }

    @Override
    public Collection<Map<String, Object>> getPointData(String chr, int pos) {
        return Collections.singleton(getFeatures(chr, pos));
    }

    public Map<String, Object> getFeatures(String chr, int start) {
        Map<String, Object> output = new HashMap<>();

        Iterator it = readFeatures(chr, start, start + 1);
        String n;

        if (it != null) {
            try {
                n = it.next();
                while (n != null) {
                    String[] scanner = n.split("\t");

                    // String chrString = scanner[0];
                    // String gff2 = scanner.next();
                    int gffStart = Integer.valueOf(scanner[3]);
                    // String gffEnd = scanner.nextInt();
                    // String gff5 = scanner.next();
                    // String gffStrand = scanner.next();
                    // String gff7 = scanner.next();
                    String gffFeature = scanner[8];
                    // scanner.close();

                    Map<String, String> valueMap = mapSplitter.split(gffFeature);

                    if (columnNames != null && !columnNames.isEmpty()) {
                        for (String name : columnNames) {
                            if (valueMap.containsKey(name)) {
                                String value = valueMap.get(name);
                                value = value.trim();
                                value = value.replaceAll("[\\s]+", "_");

                                if (output.containsKey(name)) {
                                    output.put(name, output.get(name) + "|" + value);
                                } else {
                                    output.put(name, value);
                                }
                            }
                        }
                    } else {
                        for (String name : valueMap.keySet()) {
                            String value = valueMap.get(name);
                            value = value.trim();
                            value = value.replaceAll("[\\s]+", "_");

                            if (output.containsKey(name)) {

                                output.put(name, output.get(name) + "|" + value);
                            } else {
                                output.put(name, value);
                            }
                        }
                    }
                    output.put("position", gffStart);
                    n = it.next();
                }

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        } else {
            log.trace("Strange");
        }
        // TODO Auto-generated method stub
        return output;
    }

    public Multimap<Integer, Map<String, Object>> getFeatures(String chr, int start, int end) {
        Multimap<Integer, Map<String, Object>> output = ArrayListMultimap.create();

        Iterator it = readFeatures(chr, start, end);
        String n;

        if (it != null) {
            try {
                n = it.next();
                while (n != null) {
                    String[] scanner = n.split("\t");

                    // String chrString = scanner[0];
                    // String gff2 = scanner.next();
                    int gffStart = Integer.valueOf(scanner[3]);
                    // String gffEnd = scanner.nextInt();
                    // String gff5 = scanner.next();
                    // String gffStrand = scanner.next();
                    // String gff7 = scanner.next();
                    String gffFeature = scanner[8];
                    // scanner.close();

                    Map<String, String> valueMap = mapSplitter.split(gffFeature);
                    Map<String, Object> filteredMap = new HashMap<>();

                    for (String name : columnNames) {
                        if (valueMap.containsKey(name)) {
                            if (filteredMap.containsKey(name)) {
                                filteredMap.put(name, filteredMap.get(name) + "|" + valueMap.get(name));
                            } else {
                                filteredMap.put(name, valueMap.get(name));
                            }

                        }

                    }
                    output.put(gffStart, filteredMap);
                    n = it.next();
                }

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        } else {
            log.trace("Strange");
        }
        // TODO Auto-generated method stub
        return output;
    }

    @Override
    public void setSelectAttributes(Collection<String> attributes) {
        columnNames.clear();
        columnNames.addAll(attributes);

    }

    @Override
    public Collection<VCFInfoHeaderLine> getInfoHeaderLines() {
        return Collections.emptySet();
    }

    @Override
    public VariantContext annotate(VariantContext vc) {
        VariantContextBuilder vcb = new VariantContextBuilder(vc);
        Map<String, Object> attributes = vc.getAttributes();

        Map<String, Object> features = getFeatures(vc.getChr(), vc.getStart());

        for (Entry<String, Object> entry : features.entrySet()) {
            String key = entry.getKey();

            switch (key) {
            default:
                String newKey = infoPrefix + key + infoSuffix;

                if (!attributes.containsKey(newKey)) {
                    vcb.attribute(newKey, entry.getValue());
                } else {
                    log.debug("Key is duplicated: " + newKey);
                }

                break;
            }
        }

        return vcb.make();
    }

    @Override
    public VcfLine annotate(VcfLine vc) {

        Map<String, Object> features = getFeatures(vc.getChr(), vc.getPosition());

        for (Entry<String, Object> entry : features.entrySet()) {
            String key = entry.getKey();

            switch (key) {
            default:
                String newKey = infoPrefix + key + infoSuffix;

                if (!vc.getInfo().containsKey(newKey)) {
                    vc.getInfo().put(newKey, entry.getValue());
                } else {
                    log.debug("Key is duplicated: " + newKey);
                }

                break;
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

    private String infoSuffix;

    @Override
    public void setSuffix(String infoSuffix) {
        this.infoSuffix = infoSuffix;

    }
}
