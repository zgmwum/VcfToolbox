package com.cloudinside.bio.VcfToolbox.helpers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

/**
 * Can split EFF[*].EFFECT EFF[*].IMPACT EFF[*].GENE EFF[*].CODING EFF[*].AA
 * EFF[*].TRID LOF[*].GENE NM D[*].GENE to multiple rows
 * 
 * @author pstawinski
 * 
 */
public class SplitEffTranscriptRows {

    private static final String EFF_AA_NEW_NAME = "EFF_AA";
    private static final String EFF_TRID = "EFF____TRID";
    private static final String EFF_AA = "EFF____AA";
    private static final String EFF_EFFECT = "EFF____EFFECT";
    private static final String EFF_IMPACT = "EFF____IMPACT";
    private static final String EFF_GENE = "EFF____GENE";
    private static final String EFF_CODING = "EFF____CODING";

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(SplitEffTranscriptRows.class);
    protected static final Splitter tabSplitter = Splitter.on('\t').trimResults();
    protected static final Splitter commaSplitter = Splitter.on(',').trimResults();
    protected static final Joiner tabJoiner = Joiner.on('\t');
    protected static final Joiner commaJoiner = Joiner.on(',');

    @Parameter(names = "--input", description = "Tsv input file ( - for stdin)", required = true)
    private String inputFile;
    @Parameter(names = "--output", description = "Tsv output ( - for stdout)", required = true)
    private String outputFile;

    @Parameter(names = "--effect", description = "Split by effect and gene", required = false)
    private boolean effect = false;
    @Parameter(names = "--gene", description = "Split by gene only", required = false)
    private boolean gene = false;

    // @Parameter(names = "--joinEffAndImpact", description =
    // "Join effect and impact fields", required = true)
    // private String outputFile;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        SplitEffTranscriptRows app = new SplitEffTranscriptRows();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private void go() throws IOException {
        BufferedReader br;
        if (inputFile.equals("-")) {
            br = new BufferedReader(new InputStreamReader(System.in));
        } else {
            br = new BufferedReader(new FileReader(inputFile));
        }

        BufferedWriter bw;
        if (outputFile.equals("-")) {
            bw = new BufferedWriter(new OutputStreamWriter(System.out));
        } else {
            bw = new BufferedWriter(new FileWriter(outputFile));
        }

        String line;
        List<String> header = new ArrayList<>(tabSplitter.splitToList(br.readLine()));
        Map<String, Integer> oldHeaderEntryToIndex = new HashMap<>();
        {
            int headerIndex = 0;
            for (String headerEntry : header) {
                oldHeaderEntryToIndex.put(headerEntry, headerIndex);
                headerIndex++;
            }
        }

        List<String> newHeader = new ArrayList<>(header);
        newHeader.remove(EFF_AA);
        newHeader.remove(EFF_TRID);
        int newIndex = Math.min(oldHeaderEntryToIndex.get(EFF_AA), oldHeaderEntryToIndex.get(EFF_TRID));
        newHeader.add(newIndex, EFF_AA_NEW_NAME);

        bw.write(tabJoiner.join(newHeader) + "\n");

        while ((line = br.readLine()) != null) {
            List<String> lineSplitted = new ArrayList<>(tabSplitter.splitToList(line));
            VcfLineEntry snpEffEntry = new VcfLineEntry(lineSplitted, header);
            for (List<String> outputLineList : snpEffEntry.getSplitted(newHeader)) {
                bw.write(tabJoiner.join(outputLineList) + "\n");
            }
        }
        br.close();
        bw.close();
    }

    private static class VcfLineEntry {

        private final LinkedHashMap<String, String> values = new LinkedHashMap<>();

        public VcfLineEntry(List<String> lineSplitted, List<String> header) {
            for (int i = 0; i < header.size(); i++) {
                values.put(header.get(i), lineSplitted.get(i));
            }
        }

        public List<List<String>> getSplitted(List<String> newHeader) {
            List<List<String>> result = new ArrayList<>();
            Map<String, GenePointDescription> geneToDescription = new HashMap<>();

            List<String> effects = commaSplitter.splitToList(values.get(EFF_EFFECT));
            List<String> impacts = commaSplitter.splitToList(values.get(EFF_IMPACT));
            List<String> genes = commaSplitter.splitToList(values.get(EFF_GENE));
            List<String> codings = commaSplitter.splitToList(values.get(EFF_CODING));
            List<String> aas = commaSplitter.splitToList(values.get(EFF_AA));
            List<String> trids = commaSplitter.splitToList(values.get(EFF_TRID));
            // List<String> lofs =
            // commaSplitter.splitToList(values.get("LOF[*].GENE"));
            // List<String> nmds =
            // commaSplitter.splitToList(values.get("NMD[*].GENE"));

            for (int i = 0; i < genes.size(); i++) {
                String gene = genes.get(i);
                GenePointDescription description = geneToDescription.get(gene);
                if (description == null) {
                    description = new GenePointDescription();
                    geneToDescription.put(gene, description);
                }
                try {
                    description.addLine(effects.get(i), impacts.get(i), codings.get(i), aas.get(i), trids.get(i));
                } catch (Exception e) {
                    System.out.println("here");
                }
                // ,
                // lofs.get(i), nmds.get(i)

            }

            for (Entry<String, GenePointDescription> entry : geneToDescription.entrySet()) {
                String gene = entry.getKey();
                GenePointDescription description = entry.getValue();

                LinkedHashMap<String, String> localValues = new LinkedHashMap<>(values);
                localValues.put(EFF_EFFECT,
                        description.effectSet.isEmpty() ? "" : "," + commaJoiner.join(description.effectSet) + ",");
                localValues.put(EFF_IMPACT, description.getImpact());
                localValues.put(EFF_GENE, gene);
                localValues.put(EFF_CODING, commaJoiner.join(description.codingSet));
                localValues.put(EFF_AA_NEW_NAME, description.getAA());
                // localValues.put("LOF[*].GENE",
                // commaJoiner.join(description.lofSet));
                // localValues.put("NMD[*].GENE",
                // commaJoiner.join(description.nmdSet));

                List<String> geneResult = new ArrayList<>();
                for (String header : newHeader) {
                    geneResult.add(localValues.get(header));
                }
                result.add(geneResult);

            }
            return result;
        }
    }

    private static class GenePointDescription {

        public GenePointDescription() {
        }

        public void addLine(String effect, String impact, String coding, String aa, String trid// ,
                                                                                               // String
                                                                                               // lof,
                                                                                               // String
                                                                                               // nmd
        ) {
            effectSet.add(effect);
            impactSet.add(impact);
            codingSet.add(coding);
            aaAndTridSet.add(trid + ":" + aa);
            // lofSet.add(lof);
            // nmdSet.add(nmd);
        }

        public String getImpact() {
            if (impactSet.contains("HIGH"))
                return "HIGH";
            else if (impactSet.contains("MODERATE"))
                return "MODERATE";
            else if (impactSet.contains("MODIFIER"))
                return "MODIFIER";
            else if (impactSet.contains("LOW"))
                return "LOW";
            else
                return commaJoiner.join(impactSet);
        }

        public String getAA() {
            if (aaAndTridSet.size() == 1 && aaAndTridSet.contains(".:."))
                return "";
            else
                return commaJoiner.join(aaAndTridSet);
        }

        private Set<String> effectSet = new LinkedHashSet<>();
        private Set<String> impactSet = new LinkedHashSet<>();
        private Set<String> codingSet = new LinkedHashSet<>();
        private Set<String> aaAndTridSet = new LinkedHashSet<>();
        // private Set<String> lofSet = new LinkedHashSet<>();
        // private Set<String> nmdSet = new LinkedHashSet<>();
    }
}
