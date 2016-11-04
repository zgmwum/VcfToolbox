package com.cloudinside.bio.VcfToolbox;

import java.io.Closeable;
import java.util.Collection;
import java.util.Map;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.cloudinside.bio.model.vcf.VcfLine;

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
