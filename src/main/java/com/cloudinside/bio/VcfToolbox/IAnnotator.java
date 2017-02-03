package com.cloudinside.bio.VcfToolbox;

import java.io.Closeable;
import java.util.Collection;
import java.util.Map;

import com.cloudinside.bio.model.vcf.VcfLine;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public interface IAnnotator extends Closeable {
    void setSelectAttributes(Collection<String> attributes);

    /**
     * In case of multiple annotations joins them using ','
     * 
     * @param vc
     * @return
     */
    VariantContext annotate(VariantContext vc);

    VcfLine annotate(VcfLine vc);

    Collection<VariantContext> annotateMultiplePossible(Collection<VariantContext> vcs);

    /**
     * prefix added to each column
     * 
     * @param infoPrefix
     */
    void setPrefix(String infoPrefix);

    Collection<VCFInfoHeaderLine> getInfoHeaderLines();

    Collection<Map<String, Object>> getPointData(String chr, int pos);

    void setSuffix(String infoSuffix);

}
