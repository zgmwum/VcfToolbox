package com.cloudinside.bio.VcfToolbox;

import java.util.Comparator;
import java.util.stream.Collectors;

import com.google.common.collect.ComparisonChain;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;

public class VariantContextComparatorAltAware implements Comparator<VariantContext> {
    private final VariantContextComparator vcc;

    public VariantContextComparatorAltAware(VariantContextComparator vcc) {
        this.vcc = vcc;
    }

    @Override
    public int compare(VariantContext o1, VariantContext o2) {
        return ComparisonChain.start().compare(o1, o2, vcc)
                .compare(o1.getReference().getBaseString(), o2.getReference().getBaseString())
                .compare(o1.getAlleles().stream().map(Allele::getBaseString).collect(Collectors.joining(",")),
                        o2.getAlleles().stream().map(Allele::getBaseString).collect(Collectors.joining(",")))
                .result();
    }

}
