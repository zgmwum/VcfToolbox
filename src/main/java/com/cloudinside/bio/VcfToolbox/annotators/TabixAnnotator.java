package com.cloudinside.bio.VcfToolbox.annotators;

import java.io.Closeable;
import java.io.IOException;

import org.broad.tribble.readers.TabixReader;
import org.broad.tribble.readers.TabixReader.Iterator;

public abstract class TabixAnnotator implements Closeable {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(TabixAnnotator.class);

    protected final TabixReader tabixReader;
    protected final boolean stripChrFromChromosomeNames;
    private final String filename;

    /**
     * @param filename
     *            is a tabix prepared file
     */
    public TabixAnnotator(String filename) {
        this(filename, false);

    }

    /**
     * @param filename
     *            is a tabix prepared file
     */
    public TabixAnnotator(String filename, boolean stripChrFromChromosomeNames) {
        TabixReader tabixReaderTmp = null;
        this.stripChrFromChromosomeNames = stripChrFromChromosomeNames;
        this.filename = filename;

        try {
            tabixReaderTmp = new TabixReader(filename);
        } catch (IOException e) {
            e.printStackTrace();
        }
        tabixReader = tabixReaderTmp;
    }

    public void close() throws IOException {
        tabixReader.close();
    }

    protected Iterator readFeatures(String chr, int start, int end) {
        try {
            if (stripChrFromChromosomeNames) {
                if (chr.startsWith("chr")) {
                    chr = chr.substring(3);
                }
            }

            Integer chrId = tabixReader.mChr2tid.get(chr);
            if (chrId == null) {
                log.debug("No chromosome named " + chr + " registered in " + filename);
                return null;
            } else
                return tabixReader.query(chrId, start - 1, end - 1);

            // return tabixReader.query(chr + ":" + start + "-" + end);

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

}
