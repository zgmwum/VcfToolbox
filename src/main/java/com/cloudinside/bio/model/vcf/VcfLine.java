package com.cloudinside.bio.model.vcf;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

public class VcfLine {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(VcfLine.class);

    private String chr;
    private Integer position;
    private String id;
    private String ref;
    private List<String> alt;
    private Double qual;
    private List<String> filter;
    private LinkedHashMap<String, Object> info;
    private List<String> format;
    /**
     * actually, the format encoded data; SAMPLE_NAME->FORMAT_ELEMENT->STRING
     * VALUE
     */
    private LinkedHashMap<String, Map<String, String>> samplesData;

    public VcfLine() {
        // TODO Auto-generated constructor stub
    }

    private final static Splitter commaSplitter = Splitter.on(',');
    private final static Splitter semicolonSplitter = Splitter.on(';');
    private final static Splitter colonSplitter = Splitter.on(':');

    private final static Joiner commaJoiner = Joiner.on(',');
    private final static Joiner semicolonJoiner = Joiner.on(';');
    private final static Joiner colonJoiner = Joiner.on(':');

    public static VcfLine fromLine(String line, String[] sampleNames) {
        try {
            VcfLine vcfLine = new VcfLine();
            String[] splitted = line.split("\t");
            vcfLine.setChr(splitted[0]);
            vcfLine.setPosition(Integer.valueOf(splitted[1]));
            vcfLine.setId(StringUtils.isNotBlank(splitted[2]) && !".".equals(splitted[2]) ? splitted[2] : null);
            vcfLine.setRef(splitted[3]);
            vcfLine.setAlt(commaSplitter.splitToList(splitted[4]));

            vcfLine.setQual(!".".equals(splitted[5]) ? Double.valueOf(splitted[5]) : null);
            vcfLine.setFilter(semicolonSplitter.splitToList(splitted[6]));

            LinkedHashMap<String, Object> info = new LinkedHashMap<String, Object>();
            for (String infoPart : semicolonSplitter.split(splitted[7])) {
                int equalSignIndex = infoPart.indexOf('=');
                if (equalSignIndex == -1) {
                    info.put(infoPart, "");
                } else {
                    info.put(StringUtils.substringBefore(infoPart, "="), StringUtils.substringAfter(infoPart, "="));
                }
            }

            vcfLine.setInfo(info);

            if (splitted.length > 8) {
                vcfLine.setFormat(new ArrayList<>(colonSplitter.splitToList(splitted[8])));
            }
            if (splitted.length > 9) {

                LinkedHashMap<String, Map<String, String>> samplesDataTemp = new LinkedHashMap<>();
                int index = 9;
                for (String sampleName : sampleNames) {
                    List<String> sampleData = new ArrayList<>();
                    Map<String, String> currentSampleFormat = new LinkedHashMap<>();
                    sampleData = colonSplitter.splitToList(splitted[index]);

                    int formatIndex = 0;

                    for (String formatElement : sampleData) {
                        currentSampleFormat.put(vcfLine.getFormat().get(formatIndex), formatElement);
                        formatIndex++;
                    }

                    samplesDataTemp.put(sampleName, currentSampleFormat);
                    index++;
                }

                vcfLine.setSamplesData(samplesDataTemp);
            }

            return vcfLine;
        } catch (Exception e) {
            e.printStackTrace();
            log.error("Input line: " + line);

            return null;
        }

    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public Integer getPosition() {
        return position;
    }

    public void setPosition(Integer position) {
        this.position = position;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getRef() {
        return ref;
    }

    public void setRef(String ref) {
        this.ref = ref;
    }

    public void setAlt(List<String> alt) {
        this.alt = alt;
    }

    public List<String> getAlt() {
        return alt;
    }

    public Double getQual() {
        return qual;
    }

    public void setQual(Double qual) {
        this.qual = qual;
    }

    public void setFilter(List<String> filter) {
        this.filter = filter;
    }

    public List<String> getFilter() {
        return filter;
    }

    public LinkedHashMap<String, Object> getInfo() {
        return info;
    }

    public void setInfo(LinkedHashMap<String, Object> info) {
        this.info = info;
    }

    public List<String> getFormat() {
        return format;
    }

    public void setFormat(List<String> format) {
        this.format = format;
    }

    public LinkedHashMap<String, Map<String, String>> getSamplesData() {
        return samplesData;
    }

    public void setSamplesData(LinkedHashMap<String, Map<String, String>> samplesData) {
        this.samplesData = samplesData;
    }

    public String toString(String[] samplesLine) {

        StringBuilder sb = new StringBuilder();
        try {
            sb.append(chr);
            sb.append("\t");
            sb.append(position);
            sb.append("\t");
            sb.append(dotInsteadOfNull(id));
            sb.append("\t");
            sb.append(dotInsteadOfNull(ref));
            sb.append("\t");
            sb.append(commaJoiner.join(alt));
            sb.append("\t");
            sb.append(dotInsteadOfNull(qual));
            sb.append("\t");
            sb.append(semicolonJoiner.join(filter));
            sb.append("\t");

            boolean first = true;
            for (Entry<String, Object> entry : info.entrySet()) {
                if (first) {
                    first = false;
                } else {
                    sb.append(";");
                }
                if (StringUtils.isNotBlank(entry.getValue().toString())) {
                    sb.append(entry.getKey());
                    sb.append("=");
                    sb.append(entry.getValue());
                } else {
                    sb.append(entry.getKey());
                }
            }
            if (format != null && !format.isEmpty()) {
                sb.append("\t");
                sb.append(colonJoiner.join(format));

                if (samplesLine != null) {
                    for (String sample : samplesLine) {
                        sb.append("\t");
                        sb.append(colonJoiner.join(samplesData.get(sample).values()));

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            log.error(sb.toString());
            log.error("Error in " + chr + "\t" + position);
        }

        return sb.toString();
    }

    private Object dotInsteadOfNull(Object value) {
        if (value == null) {
            return ".";
        } else {
            return value;
        }
    }

    public VcfLine(VcfLine vc) {
        super();
        this.chr = vc.chr;
        this.position = vc.position;
        this.id = vc.id;
        this.ref = vc.ref;
        this.alt = new ArrayList<>(vc.alt);
        this.qual = vc.qual;
        this.filter = new ArrayList<>(vc.filter);
        this.info = new LinkedHashMap<>(vc.info);
        if (vc.format != null)
            this.format = new ArrayList<>(vc.format);
        if (vc.samplesData != null)
            this.samplesData = new LinkedHashMap<>(vc.samplesData);
    }

}
