package com.cloudinside.bio.VcfToolbox;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.Locale;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

import org.broadinstitute.variant.vcf.VCFHeader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;
import com.google.common.collect.ImmutableSet;

/**
 * Input: single observations vcf; output: normalized observations, not merged;
 * GT updated accordingly
 * 
 * @author pstawinski
 * 
 */
public class NormalizeAllelesStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(NormalizeAllelesStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix). Cann vcf file with snpEff annotations", required = true)
    private String inputVcfFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;

    @Parameter(names = "--ref", description = "Vcf output", required = true)
    private String referenceFile;

    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    private NumberFormat formatter = new DecimalFormat("#0.00");

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        NormalizeAllelesStream app = new NormalizeAllelesStream();
        JCommander jc = null;
        try {
            jc = new JCommander(app, args);

            app.go();
        } catch (Exception e) {
            e.printStackTrace();
            jc.usage();

        }
    }

    private final static int TRANSCRIPT_NAME_POSITION = 10;
    private final static Set<String> differentFieldsThatCanBeOmmited = ImmutableSet.copyOf(new String[] { "AF", "AC" });

    private void go() throws IOException {

        ReferenceSequenceFile referenceReader = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(
                referenceFile), true);

        // TODO Auto-generated method stub
        InputStream inputStream;
        if ("-".equals(inputVcfFile)) {
            inputStream = System.in;
        } else {
            inputStream = new FileInputStream(inputVcfFile);
        }

        if (ungzip || inputVcfFile.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
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

        vcfStreamWriter.setVcfHeader(header);
        vcfStreamWriter.setSamplesLine(vcfStreamReader.getSamplesLine());

        VcfStreamReader it = vcfStreamReader;

        {
            VcfLine vc;

            while ((vc = it.next()) != null) {
                String ref = vc.getRef();
                String alt = vc.getAlt().get(0);
                int position = vc.getPosition();

                boolean changed = true;
                while (changed) {
                    changed = false;
                    if (ref.charAt(ref.length() - 1) == alt.charAt(alt.length() - 1)) {
                        // truncate rightmost if equal
                        ref = ref.substring(0, ref.length() - 1);
                        alt = alt.substring(0, alt.length() - 1);
                        changed = true;
                    }
                    if (ref.isEmpty() || alt.isEmpty()) {
                        // add base to the left, if become empty
                        position--;
                        String newNucl = new String(referenceReader.getSubsequenceAt(vc.getChr(), position, position)
                                .getBases()).toUpperCase();
                        ref = newNucl + ref;
                        alt = newNucl + alt;
                        changed = true;
                    }
                }

                while ((ref.charAt(0) == alt.charAt(0)) && (ref.length() >= 2) && (alt.length() >= 2)) {
                    ref = ref.substring(1);
                    alt = alt.substring(1);
                    position++;
                }

                vc.setPosition(position);
                vc.setRef(ref);
                vc.setAlt(Collections.singletonList(alt));
                vcfStreamWriter.write(vc);
            }
        }

        vcfStreamReader.close();
        vcfStreamWriter.close();

    }

}
