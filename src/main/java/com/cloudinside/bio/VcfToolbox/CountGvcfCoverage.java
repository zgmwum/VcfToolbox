package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VcfToolbox.utils.ChromosomePosition;
import com.google.common.base.Splitter;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class CountGvcfCoverage {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(CountGvcfCoverage.class);

    @Parameter(names = "--input", description = "GVcf input file (can be bgzipped), can be - for stdin", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Output file for statistics, may be - for stdout", required = true)
    private String output;
    @Parameter(names = "--bed", description = "Bed file with target", required = true)
    private String bedFile;

    private NumberFormat nf = new DecimalFormat("#0.0000");

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        CountGvcfCoverage app = new CountGvcfCoverage();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private final RangeSet<ChromosomePosition> target = TreeRangeSet.create();
    private final Multiset<Integer> gqMultiset = HashMultiset.create();
    private long targetSize = 0;

    private void go() throws IOException {

        {
            InputStream annotationsStream;
            if (bedFile.endsWith(".gz")) {
                annotationsStream = new GZIPInputStream(new FileInputStream(bedFile));
            } else {
                annotationsStream = new FileInputStream(bedFile);
            }
            BufferedReader br = new BufferedReader(new InputStreamReader(annotationsStream));
            String line;
            while ((line = br.readLine()) != null) {
                if (StringUtils.isNotBlank(line) && !line.startsWith("#")) {
                    String[] splitted = line.split("\t");

                    String chr = splitted[0];
                    if (chr.startsWith("chr")) {
                        chr = StringUtils.removeStart(chr, "chr");
                    }
                    chr = chr.intern();

                    Integer start = Integer.valueOf(splitted[1]);
                    Integer end = Integer.valueOf(splitted[2]);
                    targetSize += (end - start);

                    target.add(Range.openClosed(new ChromosomePosition(chr, start), new ChromosomePosition(chr, end)));
                }
            }
            IOUtils.closeQuietly(br);
        }

        InputStream inputStream;
        if ("-".equals(inputVcfFile)) {
            inputStream = System.in;
        } else {
            inputStream = new FileInputStream(inputVcfFile);
        }

        if (inputVcfFile.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
        }
        BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));

        OutputStream outputStream;
        if ("-".equals(output)) {
            outputStream = System.out;
        } else {
            outputStream = new FileOutputStream(output);
        }
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(outputStream));

        String line;

        Pattern endPattern = Pattern.compile("END=([0-9][0-9]*)");
        Splitter colonSplitter = Splitter.on(':').trimResults();
        int lineNumber = 0;

        while ((line = br.readLine()) != null) {
            lineNumber++;
            if ((lineNumber % 100000) == 0)
                System.err.print(lineNumber + " ");
            if (StringUtils.isNotBlank(line) && !line.startsWith("#")) {
                String[] splitted = line.split("\t");
                String chr = splitted[0];
                if (chr.startsWith("chr")) {
                    chr = StringUtils.removeStart(chr, "chr");
                }
                Integer start = Integer.valueOf(splitted[1]);
                String infoFields = splitted[7];
                Integer end = null;
                Matcher matcher = endPattern.matcher(infoFields);
                if (matcher.find()) {
                    end = Integer.valueOf(matcher.group(1));
                }
                String formatDeclaration = splitted[8];
                String format = splitted[9];
                int GQ_index = colonSplitter.splitToList(formatDeclaration).indexOf("GQ");
                Integer GQ_value = null;
                if (GQ_index == -1) {
                    log.warn("No GQ in " + line);
                } else {
                    try {
                        GQ_value = Integer.valueOf(colonSplitter.splitToList(format).get(GQ_index));
                    } catch (Exception e) {
                        log.error("Error parsing: " + format + "to get GQ with index: " + GQ_index);
                        e.printStackTrace();
                    }
                }

                if (GQ_value == null) {
                    log.warn("No GQvalue in " + line);
                } else if (GQ_value > 0) {
                    if (end != null) {
                        for (Integer pos = start; pos <= end; pos++) {
                            ChromosomePosition chromosomePosition = new ChromosomePosition(chr, pos);
                            register(GQ_value, chromosomePosition);
                        }
                    } else {
                        ChromosomePosition chromosomePosition = new ChromosomePosition(chr, start);
                        register(GQ_value, chromosomePosition);
                    }
                }

            }
        }

        System.err.println();

        bw.write("Target size:\t" + targetSize + "\n");
        bw.write("ge10:\t" + nf.format(greaterOrEqual(10)) + "\n");
        bw.write("ge20:\t" + nf.format(greaterOrEqual(20)) + "\n");
        bw.write("ge30:\t" + nf.format(greaterOrEqual(30)) + "\n");
        bw.write("ge40:\t" + nf.format(greaterOrEqual(40)) + "\n");
        bw.write("ge50:\t" + nf.format(greaterOrEqual(50)) + "\n");

        IOUtils.closeQuietly(bw);
        IOUtils.closeQuietly(br);

    }

    private void register(Integer GQ_value, ChromosomePosition chromosomePosition) {
        if (target.contains(chromosomePosition)) {
            if (GQ_value <= 100)
                gqMultiset.add(GQ_value);
            else
                gqMultiset.add(100);
        }
    }

    private double greaterOrEqual(int threshold) {
        int sum = 0;
        for (int i = threshold; i <= 100; i++) {
            sum += gqMultiset.count(i);
        }
        return 100. * (double) sum / (double) targetSize;
    }
}
