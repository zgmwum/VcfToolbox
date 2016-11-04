package com.cloudinside.bio;

import org.apache.commons.lang3.StringUtils;

import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;

public class VariantHashCounter {
    static private final HashFunction hf = Hashing.md5();

    /**
     * @param chr
     * @param start
     * @param end
     *            - ommited
     * @param ref
     * @param obs
     * @return
     */
    static public String hash(String chr, Integer start, Integer end, String ref, String obs) {
        StringBuilder sb = new StringBuilder();
        sb.append(chr).append(':').append(StringUtils.leftPad(Integer.toString(start), 9, '0')).append('-').append(ref)
                .append('>').append(obs);
        if (sb.length() > 100) {
            return chr + ":" + start + "-" + hf.hashUnencodedChars(sb.toString()).toString();
        } else {
            return sb.toString();
        }

    }
}
