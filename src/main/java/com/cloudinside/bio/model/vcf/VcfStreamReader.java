package com.cloudinside.bio.model.vcf;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class VcfStreamReader implements Closeable {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(VcfStreamReader.class);

    private BufferedReader bufferedReader;
    private String[] samplesLine;

    private VCFHeader vcfHeader;

    public VcfStreamReader(BufferedReader bufferedReader) {
        this.bufferedReader = bufferedReader;
    }

    private final static Pattern HEADER_LINE_PATTERN = Pattern.compile("##([^=]*)=(.*)");

    public VCFHeader getVcfHeader() {
        if (vcfHeader == null) {
            if (bufferedReader != null) {
                Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
                String line;
                try {
                    while ((line = bufferedReader.readLine()).startsWith("##")) {
                        Matcher m = HEADER_LINE_PATTERN.matcher(line);
                        if (m.matches()) {
                            String key = m.group(1);
                            String value = m.group(2);
                            headerLines.add(new VCFHeaderLine(key, value));
                        } else {
                            log.warn("Header line not matching pattern: " + line);
                        }
                    }
                    samplesLine = StringUtils.substringAfter(line, "FORMAT").trim().split("\t");
                    if (samplesLine.length == 1 && StringUtils.isBlank(samplesLine[0])) {
                        samplesLine = null;
                    }

                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
                vcfHeader = new VCFHeader(headerLines);

            }

        }
        return vcfHeader;
    }

    public String[] getSamplesLine() {
        return samplesLine;
    }

    public VcfLine next() {
        if (bufferedReader != null) {
            if (samplesLine == null) {
                getVcfHeader();
            }

            try {
                String line = bufferedReader.readLine();
                if (line == null) {
                    close();
                    return null;
                }
                VcfLine vcfLine = VcfLine.fromLine(line, samplesLine);
                return vcfLine;

            } catch (IOException e) {

                e.printStackTrace();
            }

        }
        return null;
    }

    @Override
    public void close() throws IOException {
        IOUtils.closeQuietly(bufferedReader);
        bufferedReader = null;

    }

    public static void main(String[] args) {
        try {
            VcfStreamReader x = new VcfStreamReader(new BufferedReader(
                    new FileReader("/wum/pio/experiments/wes/Sample_1ROGOWIEC/Sample_1ROGOWIEC-hcp.va.filtered.vcf")));

            System.out.println(x.getVcfHeader());
            VcfLine line;
            while (((line = x.next()) != null)) {
                System.out.println(line);
            }

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
}
