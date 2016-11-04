package com.cloudinside.bio.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;

public class RangeMultimap<K> {
    private TreeRangeMap<Integer, List<K>> map = TreeRangeMap.create();

    public static void main(String[] args) {

        RangeMultimap<String> rangeMultimap = new RangeMultimap<String>();

        rangeMultimap.put(Range.closedOpen(0, 45), "0-45");
        rangeMultimap.put(Range.closedOpen(1, 5), "1-5");
        rangeMultimap.put(Range.closedOpen(3, 15), "3-15");
        rangeMultimap.put(Range.closedOpen(15, 20), "15-20");

        for (int i = -10; i < 50; i++) {
            System.out.println(i + "\t" + rangeMultimap.get(i));
        }

        for (int i = -10; i < 50; i++) {
            System.out.println(i + "\t" + rangeMultimap.get(i, 3));
        }
    }

    public Set<K> get(int key, int radius) {
        if (radius == 1) {
            return get(key);
        }

        RangeMap<Integer, List<K>> submap = map.subRangeMap(Range.closed(key - radius, key + radius));

        Set<K> output = new HashSet<K>();
        for (List<K> values : submap.asMapOfRanges().values()) {
            output.addAll(values);
        }

        return output;
    }

    public Set<K> get(int key) {

        List<K> v = map.get(key);
        if (v == null)
            return Collections.EMPTY_SET;

        return new HashSet<K>(v);
    }

    public void put(Range<Integer> range, K value) {
        RangeMap<Integer, List<K>> subRangeMap = map.subRangeMap(range);
        Map<Range<Integer>, List<K>> submap = subRangeMap.asMapOfRanges();
        if (submap.isEmpty()) {
            ArrayList<K> list = new ArrayList<K>();
            list.add(value);
            map.put(range, list);
        } else {

            RangeSet<Integer> notCoveredSpan = TreeRangeSet.create();
            notCoveredSpan.add(range);

            List<Pair<Range<Integer>, List<K>>> newElements = new ArrayList<Pair<Range<Integer>, List<K>>>();
            for (Entry<Range<Integer>, List<K>> mapEntry : submap.entrySet()) {
                List<K> newValue = new ArrayList<K>(mapEntry.getValue());
                newValue.add(value);

                newElements.add(new ImmutablePair<Range<Integer>, List<K>>(mapEntry.getKey(), newValue));
                notCoveredSpan.remove(mapEntry.getKey());
            }

            for (Pair<Range<Integer>, List<K>> el : newElements) {
                map.put(el.getLeft(), el.getRight());
            }

            if (!notCoveredSpan.isEmpty()) {
                for (Range<Integer> notYetCoveredSpan : notCoveredSpan.asRanges()) {
                    ArrayList<K> list = new ArrayList<K>();
                    list.add(value);
                    map.put(notYetCoveredSpan, list);
                }
            }

        }
    }

}
