package com.cloudinside.bio.VcfToolbox.utils;

import com.google.common.collect.ComparisonChain;

public class ChromosomePosition implements Comparable<ChromosomePosition> {

    String chr;
    int position;

    public ChromosomePosition() {

    }

    public ChromosomePosition(String chr, int position) {
        super();
        this.chr = chr;
        this.position = position;
    }

    public String getChr() {
        return chr;
    }

    public int getPosition() {
        return position;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((chr == null) ? 0 : chr.hashCode());
        result = prime * result + position;
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
        ChromosomePosition other = (ChromosomePosition) obj;
        if (chr == null) {
            if (other.chr != null)
                return false;
        } else if (!chr.equals(other.chr))
            return false;
        if (position != other.position)
            return false;
        return true;
    }

    @Override
    public int compareTo(ChromosomePosition that) {
        return ComparisonChain.start().compare(this.chr, that.chr).compare(this.position, that.position).result();
    }

    public static Integer distance(ChromosomePosition left, ChromosomePosition right) {
        if (left == null || right == null)
            return null;
        if (left.getChr().equals(right.getChr())) {
            return Math.abs(left.getPosition() - right.getPosition());
        } else {
            return null;
        }
    }
}
