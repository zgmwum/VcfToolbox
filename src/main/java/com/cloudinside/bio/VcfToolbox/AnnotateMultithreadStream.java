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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.cloudinside.bio.VcfToolbox.annotators.BedAnnotator;
import com.cloudinside.bio.VcfToolbox.annotators.GffFeatureAnnotator;
import com.cloudinside.bio.VcfToolbox.annotators.VcfFeatureAnnotator;
import com.cloudinside.bio.VcfToolbox.annotators.VcfRadiusAnnotator;
import com.cloudinside.bio.VcfToolbox.utils.NoSplitter;
import com.cloudinside.bio.model.vcf.VcfLine;
import com.cloudinside.bio.model.vcf.VcfStreamReader;
import com.cloudinside.bio.model.vcf.VcfStreamWriter;
import com.google.common.base.CharMatcher;
import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * 
 * @author pstawinski
 * 
 */
public class AnnotateMultithreadStream {

    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
            .getLogger(AnnotateMultithreadStream.class);

    @Parameter(names = "--input", description = "Vcf input file (can be bgzipped, may need to indexed through tabix)", required = true)
    private String inputVcfFile;
    @Parameter(names = "--annotationsFile", description = "Annotations file, must be tabix encoded; at the moment only vcf is supported (you can also check gff and bed). You can separate additional attributes with : sign; for example file.vcf:infoPrefix=1000G:noChrPrefix=true:onlyColumns=column1,column2,column3", required = true, splitter = NoSplitter.class)
    private List<String> annotationsFile;
    @Parameter(names = "--output", description = "Vcf output", required = true)
    private String outputVcfFile;
    @Parameter(names = "--ignore-missing-header", description = "Ingore missing header lines", required = false)
    private boolean ignoreMissingHeader = false;
    // @Parameter(names = "--ignore-input-index", description =
    // "Require index of input", required = false)
    private boolean ignoreInputIndex = true;
    @Parameter(names = "--threads", description = "Ingore missing header lines", required = false)
    private int threads = 4;
    @Parameter(names = "--ungzip", description = "Ungzip input", required = false)
    private boolean ungzip = false;

    private VcfStreamReader it;

    private Stopwatch stopwatch;

    private int counter;

    private VcfStreamWriter vcfWriter;

    public static void main(String[] args) {
        // to have Double formatted correctly
        Locale.setDefault(Locale.ENGLISH);

        AnnotateMultithreadStream app = new AnnotateMultithreadStream();
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

        VcfStreamReader vcfStreamReader = new VcfStreamReader(new BufferedReader(new InputStreamReader(inputStream)));

        VCFHeader header = vcfStreamReader.getVcfHeader();
        it = vcfStreamReader;

        OutputStream outputStream;
        if ("-".equals(outputVcfFile)) {
            outputStream = System.out;
        } else {
            outputStream = new FileOutputStream(outputVcfFile);
        }
        vcfWriter = new VcfStreamWriter(new BufferedWriter(new OutputStreamWriter(outputStream)));

        List<IAnnotator> annotators = openAnnotators(annotationsFile);

        // build header

        for (IAnnotator annotator : annotators) {
            for (VCFInfoHeaderLine line : annotator.getInfoHeaderLines()) {
                if (header.getInfoHeaderLine(line.getID()) == null) {
                    header.addMetaDataLine(line);
                } else {
                    log.debug("Ommiting line, as already there " + line);
                }
            }
        }

        for (IAnnotator annotator : annotators) {
            IOUtils.closeQuietly(annotator);
        }

        vcfWriter.setSamplesLine(vcfStreamReader.getSamplesLine());
        vcfWriter.setVcfHeader(header);

        counter = 0;
        stopwatch = Stopwatch.createStarted();

        ExecutorService executor = Executors.newFixedThreadPool(threads);
        for (int i = 0; i < threads; i++) {
            Worker worker = new Worker(this, annotationsFile);
            executor.execute(worker);
        }

        executor.shutdown();
        try {
            executor.awaitTermination(10, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        IOUtils.closeQuietly(vcfStreamReader);
        for (IAnnotator annotator : annotators) {
            IOUtils.closeQuietly(annotator);
        }
        vcfWriter.close();

    }

    protected List<IAnnotator> openAnnotators(List<String> files) {
        List<IAnnotator> annotators = new ArrayList<IAnnotator>(files.size());
        for (String annotationFile : files) {
            int delimiter = annotationFile.indexOf(':');
            String fileName;
            List<String> columns = null;
            boolean stripChrFromChromosomeName = false;
            String infoPrefix = "";
            String infoSuffix = "";

            int radius = 1;

            if (delimiter == -1) {
                fileName = annotationFile;
            } else {
                fileName = annotationFile.substring(0, delimiter);

                Map<String, String> argumentsMap = Splitter.on(":").withKeyValueSeparator("=")
                        .split(annotationFile.substring(delimiter + 1));

                for (Entry<String, String> data : argumentsMap.entrySet()) {
                    switch (data.getKey()) {
                    case "noChrPrefix":
                        if ("true".equals(data.getValue())) {
                            stripChrFromChromosomeName = true;
                        }
                        break;
                    case "infoPrefix":
                        infoPrefix = data.getValue();
                        break;
                    case "infoSuffix":
                        infoSuffix = data.getValue();
                        break;
                    case "onlyColumns":
                        columns = Arrays.asList(data.getValue().split(","));
                        break;
                    case "radius":
                        radius = Integer.valueOf(data.getValue());
                        break;
                    default:
                        log.warn("Not recognized data " + data);
                        break;
                    }

                }
            }

            IAnnotator annotator = null;
            if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
                if (radius == 1)
                    annotator = new VcfFeatureAnnotator(fileName, columns, false, stripChrFromChromosomeName);
                else
                    annotator = new VcfRadiusAnnotator(fileName, stripChrFromChromosomeName, radius);
            } else if (fileName.endsWith(".gff") || fileName.endsWith(".gff.gz")) {
                annotator = new GffFeatureAnnotator(fileName, columns);
            } else if (fileName.endsWith(".bed") || fileName.endsWith(".bed.gz")) {
                annotator = new BedAnnotator(fileName, columns.get(0),
                        CharMatcher.JAVA_LETTER_OR_DIGIT.negate().replaceFrom(fileName, "_"));
            } else {
                throw new NullPointerException("unable to detect file type from filename: " + fileName);
            }

            annotator.setPrefix(infoPrefix);
            annotator.setSuffix(infoSuffix);

            if (annotator != null) {
                annotators.add(annotator);
            }
        }
        return annotators;
    }

    private synchronized VcfLine next() {
        if (++counter % 10000 == 0) {
            long seconds = stopwatch.elapsed(TimeUnit.SECONDS);
            log.debug("Processed " + counter + ", " + (counter / seconds) + " entries / s");
        }

        return it.next();

    }

    private static class Worker implements Runnable {
        private final AnnotateMultithreadStream parent;
        private final List<IAnnotator> annotators;

        public Worker(AnnotateMultithreadStream parent, List<String> annotationsFile) {
            this.parent = parent;

            this.annotators = parent.openAnnotators(annotationsFile);

        }

        @Override
        public void run() {
            VcfLine vc;
            while ((vc = parent.next()) != null) {
                for (IAnnotator annotator : annotators) {
                    vc = annotator.annotate(vc);
                }
                parent.write(vc);
            }

            for (IAnnotator annotator : annotators) {
                IOUtils.closeQuietly(annotator);
            }
        }
    }

    public synchronized void write(VcfLine vc) {
        vcfWriter.write(vc);

    }

}
