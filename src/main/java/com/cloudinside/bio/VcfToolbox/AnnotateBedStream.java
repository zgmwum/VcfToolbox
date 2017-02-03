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
import java.util.Locale;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VcfToolbox.utils.ChromosomePosition;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;
import com.google.common.base.Optional;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Read no-sample vcf file
 * 
 * Writes file with multiallele obs splitted;
 * 
 * @author pstawinski
 * 
 */
public class AnnotateBedStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(AnnotateBedStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix), use - for stdin", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    @Parameter(names = "--annotations", description = "Bed file with annotations", required = true)
    private String bedFile;
    @Parameter(names = "--description", description = "INFO field", required = true)
    private String description;
    @Parameter(names = "--fieldName", description = "Field name", required = true)
    private String fieldName;
    @Parameter(names = "--useBedValue", description = "Field name", required = false)
    private Boolean useBedValue = Boolean.FALSE;

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        AnnotateBedStream app = new AnnotateBedStream();
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
        // TODO Auto-generated method stub

        InputStream inputStream;
        if ("-".equals(inputVcfFile)) {
            inputStream = System.in;
        } else {
            inputStream = new FileInputStream(inputVcfFile);
        }

        if (ungzip) {
            inputStream = new GZIPInputStream(inputStream);
        }

        TreeSet<ChromosomePosition> borders = new TreeSet<>();
        RangeMap<ChromosomePosition, String> annotations = TreeRangeMap.create();
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
                    String value = "";
                    Integer start = Integer.valueOf(splitted[1]);
                    Integer end = Integer.valueOf(splitted[2]);

                    if (splitted.length >= 4 && useBedValue)
                        value = splitted[3].intern();
                    else {
                        value = "0";
                    }

                    ChromosomePosition rangeStart = new ChromosomePosition(chr, start);
                    ChromosomePosition rangeEnd = new ChromosomePosition(chr, end);
                    annotations.put(Range.openClosed(rangeStart, rangeEnd), value);

                    borders.add(rangeStart);
                    borders.add(rangeEnd);
                }

            }
            IOUtils.closeQuietly(br);
        }

        VcfStreamReader vcfStreamReader = new VcfStreamReader(new BufferedReader(new InputStreamReader(inputStream)));
        VCFHeader header = vcfStreamReader.getVcfHeader();

        OutputStream outputStream;
        if ("-".equals(outputVcfFile)) {
            outputStream = System.out;
        } else {
            outputStream = new FileOutputStream(outputVcfFile);
        }
        VcfStreamWriter vcfStreamWriter = new VcfStreamWriter(new BufferedWriter(new OutputStreamWriter(outputStream)));

        header.addMetaDataLine(new VCFInfoHeaderLine(fieldName, 1, VCFHeaderLineType.String, description));
        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        int counter = 0;

        VcfLine vc;
        while ((vc = vcfStreamReader.next()) != null) {

            String chr = vc.getChr();
            if (chr.startsWith("chr")) {
                chr = StringUtils.removeStart(chr, "chr");
            }
            ChromosomePosition position = new ChromosomePosition(chr, vc.getPosition());

            if (annotations.get(position) != null) {
                vc.getInfo().put(fieldName, annotations.get(position));
            } else if (!useBedValue) {
                // use distance to the nearest annotated area
                int distanceLeft = Integer.MAX_VALUE;

                ChromosomePosition left = borders.floor(position);
                distanceLeft = Optional.fromNullable(ChromosomePosition.distance(left, position)).or(Integer.MAX_VALUE);

                int distanceRight = Integer.MAX_VALUE;

                ChromosomePosition right = borders.ceiling(position);
                distanceRight = Optional.fromNullable(ChromosomePosition.distance(right, position))
                        .or(Integer.MAX_VALUE);

                int distance = Math.min(distanceLeft, distanceRight);
                vc.getInfo().put(fieldName, distance);
            }

            vcfStreamWriter.write(vc);

            if (++counter % 10000 == 0) {
                System.err.println("Processed " + counter);
            }
        }

        IOUtils.closeQuietly(vcfStreamWriter);
        IOUtils.closeQuietly(vcfStreamReader);

    }
}
