package com.cloudinside.bio.VcfToolbox.annotators;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.broad.tribble.readers.TabixReader.Iterator;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.cloudinside.bio.VcfToolbox.IAnnotator;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.google.common.base.Splitter;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class VcfFeatureAnnotator extends TabixAnnotator implements IAnnotator {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(VcfFeatureAnnotator.class);

    private final Collection<String> columnNames;
    private final boolean populateSnpId;
    private final String filename;

    private String infoPrefix;

    /**
     * @param filename
     *            is a tabix preapred file
     */
    public VcfFeatureAnnotator(String filename, Collection<String> columnNames, boolean populateSnpId,
            boolean stripChrFromChromosomeName) {

        super(filename, stripChrFromChromosomeName);
        this.columnNames = columnNames;
        this.populateSnpId = populateSnpId;
        this.filename = filename;
    }

    @Override
    public Collection<Map<String, Object>> getPointData(String chr, int pos) {
        return getFeatures(chr, pos, pos).get(pos);
    }

    public Map<String, Object> getFeatures(String chr, int pos, String ref, String obs) {
        Map<String, Object> output = new HashMap<String, Object>();

        Iterator it = readFeatures(chr, pos, pos + 1);
        String next;

        if (it != null) {
            try {
                next = it.next();
                while (next != null) {
                    parseVcfLine(ref, obs, output, next);
                    next = it.next();
                }

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        } else {
            // log.debug("No iterator at all for " + filename + " @ " + chr +
            // " " + pos);
        }
        // TODO Auto-generated method stub
        return output;
    }

    public Multimap<Integer, Map<String, Object>> getFeatures(String chr, int start, int stop) {

        Iterator it = readFeatures(chr, start, stop);
        String next;

        Multimap<Integer, Map<String, Object>> multimap = HashMultimap.create();

        if (it != null) {
            try {
                next = it.next();
                while (next != null) {
                    Map<String, Object> output = new HashMap<String, Object>();
                    String[] line = next.split("\t");
                    int position = Integer.valueOf(line[1]);
                    String readSnpId = line[2];
                    String readRef = line[3];

                    List<String> readObs = Splitter.on(',').splitToList(line[4]);
                    String readInfo = line[7];

                    Map<String, Object> infoValues = parseMap(readInfo);

                    if (!".".equals(readSnpId))
                        output.put("ID", readSnpId);
                    for (String col : columnNames) {
                        if (infoValues.containsKey(col))
                            output.put(col, infoValues.get(col));
                    }
                    output.put("OBS", line[4]); // read obs raw
                    output.put("REF", readRef); // read obs raw

                    multimap.put(position, output);

                    next = it.next();
                }

            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        } else {
            log.debug("Strange");
        }
        // TODO Auto-generated method stub
        return multimap;
    }

    private void parseVcfLine(String ref, String obs, Map<String, Object> output, String n) {
        String[] line = n.split("\t");
        String readSnpId = line[2];
        String readRef = line[3];

        List<String> readObs = Splitter.on(',').splitToList(line[4]);
        String readInfo = line[7];

        if (ref.equalsIgnoreCase(readRef) && readObs.contains(obs)) {
            Map<String, Object> infoValues = parseMap(readInfo);

            if (!".".equals(readSnpId))
                output.put("ID", readSnpId);
            if (columnNames != null && !columnNames.isEmpty()) {
                for (String col : columnNames) {
                    output.put(col, infoValues.get(col));
                }
            } else {
                output.putAll(infoValues);
            }
            output.put("OBS", line[4]); // read obs raw
            output.put("REF", readRef); // read obs raw
        }
    }

    private Map<String, Object> parseMap(String readInfo) {
        Map<String, Object> map = new HashMap<>();
        if (StringUtils.isBlank(readInfo.trim()))
            return map;

        for (String keyValue : readInfo.split(";")) {
            String key, value = "";
            if (keyValue.indexOf('=') != -1) {
                key = keyValue.substring(0, keyValue.indexOf('='));
                value = keyValue.substring(keyValue.indexOf('=') + 1);
            } else {
                key = keyValue;
            }
            map.put(key, value);
        }

        return map;
    }

    @Override
    public void setSelectAttributes(Collection<String> attributes) {
        columnNames.clear();
        columnNames.addAll(attributes);

    }

    @Override
    public Collection<VCFInfoHeaderLine> getInfoHeaderLines() {
        VCFFileReader vcfFileReader = new VCFFileReader(new File(filename), false);
        VCFHeader header = vcfFileReader.getFileHeader();
        vcfFileReader.close();

        Set<VCFInfoHeaderLine> lines = new HashSet<>();
        for (VCFInfoHeaderLine infoHeaderLine : header.getInfoHeaderLines()) {
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

                switch (key) {
                case "ID":
                    if ((StringUtils.isBlank(vc.getID()) || ".".equals(vc.getID()))) {
                        vcb.id((String) entry.getValue());
                    }
                    break;
                case "OBS":
                case "REF":
                    // omit
                    break;
                default:

                    String newKey = infoPrefix + key + infoSuffix;
                    if (!attributes.containsKey(newKey)) {
                        vcb.attribute(newKey, entry.getValue());
                    } else {
                        log.debug("Duplicated key " + newKey);
                    }

                    break;
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

                switch (key) {
                case "ID":
                    if ((StringUtils.isBlank(vc.getId()) || ".".equals(vc.getId()))) {
                        vc.setId((String) entry.getValue());
                    }
                    break;
                case "OBS":
                case "REF":
                    // omit
                    break;
                default:

                    String newKey = infoPrefix + key + infoSuffix;
                    if (!vc.getInfo().containsKey(newKey)) {
                        vc.getInfo().put(newKey, entry.getValue());
                    } else {
                        log.debug("Duplicated key " + newKey);
                    }

                    break;
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

    private String infoSuffix;

    @Override
    public void setSuffix(String infoSuffix) {
        this.infoSuffix = infoSuffix;

    }
}
