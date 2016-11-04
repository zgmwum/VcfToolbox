package com.cloudinside.bio.model.vcf;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

public class VcfStreamWriter implements Closeable {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(VcfStreamWriter.class);

    private BufferedWriter bufferedWriter;
    private String[] samplesLine = null;

    private VCFHeader vcfHeader;
    private boolean headerWritten = false;
    private boolean hasFormatColumn;

    public VcfStreamWriter(BufferedWriter bufferedWriter) {
        this.bufferedWriter = bufferedWriter;
    }

    public void setSamplesLine(String[] samplesLine) {
        this.samplesLine = samplesLine;
        if (samplesLine != null && samplesLine.length != 0)
            hasFormatColumn = true;
    }

    public void setVcfHeader(VCFHeader vcfHeader) {
        this.vcfHeader = vcfHeader;
    }

    public void write(VcfLine vcfLine) {
        if (!headerWritten)
            writeHeader();

        try {
            bufferedWriter.write(vcfLine.toString(samplesLine));
            bufferedWriter.write("\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void writeHeader() {
        try {
            // bufferedWriter.write("##fileformat=VCFv4.1\n");
            for (VCFHeaderLine line : vcfHeader.getMetaDataInInputOrder()) {
                bufferedWriter.write("##");
                bufferedWriter.write(line.toString());
                bufferedWriter.write("\n");
            }
            bufferedWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
            if (hasFormatColumn && samplesLine != null && samplesLine.length > 0)
                bufferedWriter.write("\tFORMAT");

            if (samplesLine != null) {
                for (String sample : samplesLine) {
                    bufferedWriter.write("\t");
                    bufferedWriter.write(sample);
                }
            }
            bufferedWriter.write("\n");

        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        headerWritten = true;

    }

    @Override
    public void close() throws IOException {
        IOUtils.closeQuietly(bufferedWriter);

    }
}
