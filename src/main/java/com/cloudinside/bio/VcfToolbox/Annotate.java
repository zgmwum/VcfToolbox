package com.cloudinside.bio.VcfToolbox;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VcfToolbox.annotators.GffFeatureAnnotator;
import com.cloudinside.bio.VcfToolbox.annotators.VcfFeatureAnnotator;
import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;

/**
 * 
 * @author pstawinski
 * 
 */
public class Annotate {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(Annotate.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--annotationsFile", description = "Annotations file, must be tabix encoded; at the moment only vcf is supported", required = true)
    private List<String> annotationsFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header lines", required = false)
    private boolean ignoreMissingHeader = false;
    // @Parameter(names = "--ignore-input-index", description =
    // "Require index of input", required = false)
    private boolean ignoreInputIndex = true;
    // tabix is NOT threadsafe
    // @Parameter(names = "--threads", description =
    // "Ingore missing header lines", required = false)
    private int threads = 1;

    private CloseableIterator<VariantContext> it;

    private Stopwatch stopwatch;

    private int counter;

    private VariantContextWriter vcfWriter;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        Annotate app = new Annotate();
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
        // TODO Auto-generated method stub
        String filename = inputVcfFile; // "/archive/pio/tmp/merged_bcf.changed.vcf.gz";
        File file = new File(filename);

        VCFFileReader vcfFileReader = new VCFFileReader(file, false);
        VCFHeader header = vcfFileReader.getFileHeader();
        it = vcfFileReader.iterator();

        EnumSet<Options> options = EnumSet.noneOf(Options.class);
        if (ignoreMissingHeader)
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

        vcfWriter = VariantContextWriterFactory.sortOnTheFly(
                VariantContextWriterFactory.create(new File(outputVcfFile), header.getSequenceDictionary(), options),
                10000);

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

        // build header

        for (IAnnotator annotator : annotators) {
            for (VCFInfoHeaderLine line : annotator.getInfoHeaderLines()) {
                if (header.getInfoHeaderLine(line.getID()) == null) {
                    header.addMetaDataLine(line);
                } else {
                    log.debug("Ommiting line " + line);
                }
            }
        }

        vcfWriter.writeHeader(header);

        counter = 0;
        stopwatch = Stopwatch.createStarted();

        ExecutorService executor = Executors.newFixedThreadPool(threads);
        for (int i = 0; i < threads; i++) {
            Worker worker = new Worker(this, annotators);
            executor.execute(worker);
        }

        executor.shutdown();
        try {
            executor.awaitTermination(10, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        IOUtils.closeQuietly(vcfFileReader);
        for (IAnnotator annotator : annotators) {
            IOUtils.closeQuietly(annotator);
        }
        vcfWriter.close();

    }

    private synchronized VariantContext next() {
        if (++counter % 10000 == 0) {
            long seconds = stopwatch.elapsed(TimeUnit.SECONDS);
            log.debug("Processed " + counter + ", " + (counter / seconds) + " entries / s");
        }

        if (it.hasNext()) {
            return it.next();
        } else {
            return null;
        }
    }

    private static class Worker implements Runnable {
        private final Annotate parent;
        private final List<IAnnotator> annotators;

        public Worker(Annotate parent, List<IAnnotator> annotators) {
            this.parent = parent;
            this.annotators = annotators;

        }

        @Override
        public void run() {
            VariantContext vc;
            while ((vc = parent.next()) != null) {
                for (IAnnotator annotator : annotators) {
                    vc = annotator.annotate(vc);
                }
                parent.write(vc);
            }
        }

    }

    public synchronized void write(VariantContext vc) {
        vcfWriter.add(vc);

    }

}
