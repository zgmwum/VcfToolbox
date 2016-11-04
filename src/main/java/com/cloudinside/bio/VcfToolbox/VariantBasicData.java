package com.cloudinside.bio.VcfToolbox;

import com.cloudinside.bio.VariantHashCounter;
import com.google.common.collect.ComparisonChain;

public class VariantBasicData implements Comparable<VariantBasicData> {
    String chr, ref, obs;
    int pos;

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public String getRef() {
        return ref;
    }

    public void setRef(String ref) {
        this.ref = ref;
    }

    public String getObs() {
        return obs;
    }

    public void setObs(String obs) {
        this.obs = obs;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    public VariantBasicData() {
        // TODO Auto-generated constructor stub
    }

    public VariantBasicData(String chr, int pos, String ref, String obs) {
        super();
        this.chr = chr;
        this.ref = ref;
        this.obs = obs;
        this.pos = pos;
    }

    @Override
    public String toString() {
        return "VariantBasicData [chr=" + chr + ", ref=" + ref + ", obs=" + obs + ", pos=" + pos + "]";
    }

    @Override
    public int compareTo(VariantBasicData that) {
        return ComparisonChain.start().compare(this.chr, that.chr).compare(this.pos, that.pos)
                .compare(this.ref, that.ref).compare(this.obs, that.obs).result();

    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((chr == null) ? 0 : chr.hashCode());
        result = prime * result + ((obs == null) ? 0 : obs.hashCode());
        result = prime * result + pos;
        result = prime * result + ((ref == null) ? 0 : ref.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        VariantBasicData other = (VariantBasicData) obj;
        if (chr == null) {
            if (other.chr != null)
                return false;
        } else if (!chr.equals(other.chr))
            return false;
        if (obs == null) {
            if (other.obs != null)
                return false;
        } else if (!obs.equals(other.obs))
            return false;
        if (pos != other.pos)
            return false;
        if (ref == null) {
            if (other.ref != null)
                return false;
        } else if (!ref.equals(other.ref))
            return false;
        return true;
    }

    public String hash() {
        return VariantHashCounter.hash(chr, pos, pos, ref, obs);
    }

}
