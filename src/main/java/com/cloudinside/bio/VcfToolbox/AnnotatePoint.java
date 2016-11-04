package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VcfToolbox.annotators.GffFeatureAnnotator;
import com.cloudinside.bio.VcfToolbox.annotators.VcfFeatureAnnotator;
import com.cloudinside.bio.VcfToolbox.utils.ChromosomePosition;
import com.cloudinside.bio.VcfToolbox.utils.NoSplitter;
import com.google.common.base.Splitter;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;

/**
 * Can read file with chromosome\tposition\treference\tchange data only and
 * annotate wisely. It tries to behave correctly if there is no chromosome /
 * reference / change data.
 * 
 * @author pstawinski
 * 
 */
public class AnnotatePoint {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(AnnotatePoint.class);

    @Parameter(names = "--input", description = "Tab separated input file", required = true)
    private String inputFile;
    @Parameter(names = "--annotationsFile", description = "Annotations file, must be tabix encoded; at the moment only vcf is supported (you can also check gff and bed). You can separate additional attributes with : sign; for example file.vcf:infoPrefix=1000G:noChrPrefix=true:onlyColumns=column1,column2,column3", required = true, splitter = NoSplitter.class)
    private List<String> annotationsFile;
    @Parameter(names = "--output", description = "Output file", required = true)
    private String outputFile;

    @Parameter(names = "--chromosomes", description = "Chromosomes list to be used when there is no chromosomal data, delimited by space", required = false)
    private String chromosomes;

    /**
     * @param args
     */
    public static void main(String[] args) {
        Locale.setDefault(Locale.ENGLISH);

        AnnotatePoint app = new AnnotatePoint();
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
        try {
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));

            Splitter splitter = Splitter.on('\t').trimResults();

            List<IAnnotator> annotators = new ArrayList<IAnnotator>(annotationsFile.size());
            for (String annotationFile : annotationsFile) {
                int delimiter = annotationFile.indexOf(':');
                String fileName;
                List<String> columns;
                boolean stripChrFromChromosomeName = false;
                if (delimiter == -1) {
                    columns = null;
                    fileName = annotationFile;
                } else {
                    fileName = annotationFile.substring(0, delimiter);
                    List<String> additionalData = Splitter.on(',').splitToList(annotationFile.substring(delimiter + 1));
                    columns = new ArrayList<>();
                    for (String data : additionalData) {
                        if (data.startsWith("C-")) {
                            columns.add(data.substring(2));
                        } else if (data.equals("noprefix")) {
                            stripChrFromChromosomeName = true;
                        } else {
                            log.warn("Not recognized data " + data);
                        }
                    }
                }

                IAnnotator annotator = null;
                if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
                    annotator = new VcfFeatureAnnotator(fileName, columns, false, stripChrFromChromosomeName);
                } else if (fileName.endsWith(".gff") || fileName.endsWith(".gff.gz")) {
                    annotator = new GffFeatureAnnotator(fileName, columns);
                } else if (fileName.endsWith(".bed") || fileName.endsWith(".bed.gz")) {
                    throw new NullPointerException("bed not supported yet");
                }

                if (annotator != null) {
                    annotators.add(annotator);
                }
            }

            // EnumSet<Options> options = EnumSet.noneOf(Options.class);
            // options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

            // SAMSequenceDictionary sd = new SAMSequenceDictionary();
            // VariantContextWriter vcfWriter =
            // VariantContextWriterFactory.sortOnTheFly(
            // VariantContextWriterFactory.create(new File(outputFile), sd,
            // options), 10000);

            Map<ChromosomePosition, Multimap<String, Object>> outputMap = new LinkedHashMap<ChromosomePosition, Multimap<String, Object>>();

            Set<String> columns = new LinkedHashSet<String>();
            String line;
            while ((line = br.readLine()) != null) {
                List<String> splitted = splitter.splitToList(line);
                String chr;
                // , obs, ref;
                int position = 0;
                // if (splitted.size() > 3) {
                // obs = splitted.get(3);
                // }
                // if (splitted.size() > 2) {
                // ref = splitted.get(2);
                // }
                position = Integer.valueOf(splitted.get(1));
                chr = splitted.get(0);
                Multimap<String, Object> localMultimap = ArrayListMultimap.create();
                outputMap.put(new ChromosomePosition(chr, position), localMultimap);

                int annotatorIndex = 0;
                for (IAnnotator annotator : annotators) {
                    annotatorIndex++;
                    Collection<Map<String, Object>> annotations = annotator.getPointData(chr, position);
                    for (Map<String, Object> map : annotations) {
                        for (Entry<String, Object> entry : map.entrySet()) {
                            String key = annotatorIndex + "_" + entry.getKey();
                            localMultimap.put(key, entry.getValue());

                            columns.add(key);
                        }
                    }
                }
            }

            bw.write("chr\tpos\t" + StringUtils.join(columns, '\t') + "\n");

            for (Entry<ChromosomePosition, Multimap<String, Object>> entry : outputMap.entrySet()) {
                bw.write(entry.getKey().getChr() + "\t" + entry.getKey().getPosition());
                for (String column : columns) {
                    bw.write("\t");
                    Collection<Object> valuesCollection = entry.getValue().get(column);
                    bw.write(valuesCollection.size() > 1 ? valuesCollection.toString() : Iterables.getFirst(
                            valuesCollection, "").toString());
                }
                bw.write("\n");
            }

            IOUtils.closeQuietly(bw);
            // vcfWriter.close();
            IOUtils.closeQuietly(br);

        } catch (FileNotFoundException e) {

            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
}
