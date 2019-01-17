package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

/**
 * Translates tsv file with column headers to vcf
 * 
 * @author pstawinski f.ex.
 * 
 *         <pre>
 *  
 * --input /tmp/MutationsSomatic.txt --chr chr --pos start --ref ref --alt obs --ignore end --refseq /wum/pio/processing-data/reference-seq/rcrs/rcrs.fasta
 *         </pre>
 */
public class CsvToVcf {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(CsvToVcf.class);

    @Parameter(names = "--input", description = "Tsv file", required = true)
    private String inputTsvFile;
    @Parameter(names = "--output", description = "Tsv file", required = true)
    private String outputFile;
    @Parameter(names = "--chr", description = "Chromosome column name", required = true)
    private String chrColName;
    @Parameter(names = "--pos", description = "Position column name", required = true)
    private String posColName;
    @Parameter(names = "--ref", description = "Ref column name", required = true)
    private String refColName;
    @Parameter(names = "--alt", description = "Alt column name", required = true)
    private String altColName;
    @Parameter(names = "--refseq", description = "Reference sequence in fasta, can be used to repair missing data like: A>- instead of CA>A", required = false)
    private String refseq;
    @Parameter(names = "--commentIndicator", description = "How may the comment start", required = false)
    private List<String> commentLineStart = ImmutableList.<String> builder().add("#").add("\"#").build();

    @Parameter(names = "--ignore", description = "Column names to be ignored, add multiple times if needed", required = false)
    private List<String> ignoreColName = new ArrayList<>();

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        CsvToVcf app = new CsvToVcf();
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

        posColName = posColName.replace("\"", "").replace("-", "_");
        chrColName = chrColName.replace("\"", "").replace("-", "_");
        refColName = refColName.replace("\"", "").replace("-", "_");
        altColName = altColName.replace("\"", "").replace("-", "_");

        // TODO Auto-generated method stub
        String filename = inputTsvFile; // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        File file = new File(filename);
        BufferedReader br = null;
        BufferedWriter bw = null;
        Set<String> ignoreColNameSet = new HashSet<>(ignoreColName);
        try {
            if (filename.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
            } else {
                br = new BufferedReader(new FileReader(file));
            }

            bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputFile))));

            bw.write("##fileformat=VCFv4.1\n");

            Map<String, byte[]> refseqMap = null;
            if (refseq != null) {
                refseqMap = new HashMap<>();
                FastaSequenceFile fastaSequenceFile = new FastaSequenceFile(new File(refseq), true);
                ReferenceSequence rs;
                while ((rs = fastaSequenceFile.nextSequence()) != null) {
                    refseqMap.put(rs.getName(), rs.getBases());
                }
            }

            String line;
            boolean vcfHeaderWritten = false;
            Map<String, Integer> columnNames = null;
            List<String> columnNamesList = new ArrayList<>();
            while ((line = br.readLine()) != null) {

                boolean isComment = false;
                for (String commentInd : commentLineStart) {
                    if (line.startsWith(commentInd))
                        isComment = true;
                }

                if (isComment) {
                    // omit
                } else if (columnNames == null) {
                    columnNames = new HashMap<String, Integer>();
                    int index = 0;
                    for (String colName : Splitter.on('\t').splitToList(line)) {
                        colName = colName.trim().replace("\"", "").replace("-", "_");
                        columnNames.put(colName, index);
                        index++;
                        columnNamesList.add(colName);
                    }

                    for (String colName : columnNames.keySet()) {
                        if (ignoreColNameSet.contains(colName) || chrColName.equals(colName)
                                || posColName.equals(colName) || refColName.equals(colName)
                                || altColName.equals(colName)) {
                            // omit
                        } else {
                            bw.write("##INFO=<ID=" + colName
                                    + ",Number=.,Type=String,Description=\"Description needed\">\n");
                        }
                    }
                } else {
                    if (!vcfHeaderWritten) {
                        vcfHeaderWritten = true;
                        bw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
                    }

                    try {

                        int index = 0;
                        int position = 0;
                        String chr = null;
                        String ref = null, alt = null;
                        Map<String, String> columnToValue = new HashMap<>();
                        for (String value : Splitter.on('\t').splitToList(line)) {

                            value = value.trim().replace("\"", "").replace(" ", "_").replace(";", "_");
                            // .replace("/", "_")
                            // .replace("'", "_")
                            // .replace("&", "_");
                            String colName = columnNamesList.get(index);
                            if (ignoreColNameSet.contains(colName)) {
                                // omit
                            } else if (chrColName.equals(colName)) {
                                chr = value;
                            } else if (posColName.equals(colName)) {
                                if (".".equals(value)) {
                                    // omit
                                    log.warn("Got . position, ommiting at line " + line);
                                    position = -1;
                                } else {
                                    position = Integer.valueOf(value);
                                }
                            } else if (refColName.equals(colName)) {
                                ref = value;
                            } else if (altColName.equals(colName)) {
                                alt = value;
                            } else {
                                columnToValue.put(colName, value);
                            }
                            index++;
                        }

                        if ((StringUtils.isBlank(alt) || alt.equals("-")) && refseqMap != null) {
                            byte[] sequence = refseqMap.get(chr);

                            if (sequence == null) {
                                log.error("chromosome not recognized, variant ommited: " + chr);
                                continue;
                            }

                            String prefixBase = new String(new byte[] { sequence[position - 2] });
                            position -= 1;
                            alt = prefixBase;
                            ref = prefixBase + ref;
                        }

                        if ((StringUtils.isBlank(ref) || ref.equals("-")) && refseqMap != null) {
                            byte[] sequence = refseqMap.get(chr);

                            if (sequence == null) {
                                log.error("chromosome not recognized, variant ommited: " + chr);
                                continue;
                            }

                            String prefixBase = new String(new byte[] { sequence[position - 2] });
                            position -= 1;
                            alt = prefixBase + alt;
                            ref = prefixBase;
                        }

                        if (position != -1) {
                            bw.write(chr + "\t" + position + "\t.\t" + ref.toUpperCase() + "\t" + alt.toUpperCase()
                                    + "\t" + ".\t.\t" + Joiner.on(";").withKeyValueSeparator("=").join(columnToValue)
                                    + "\n");
                        }
                    } catch (Exception e) {
                        log.debug("Error processing line: " + line, e);
                        throw e;
                    }
                }
            }

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        IOUtils.closeQuietly(br);
        IOUtils.closeQuietly(bw);

    }
}
