package com.cloudinside.bio.VcfToolbox.utils;

import java.util.Collections;
import java.util.List;

import com.beust.jcommander.converters.IParameterSplitter;

public class NoSplitter implements IParameterSplitter {

    public List<String> split(String value) {
        return Collections.singletonList(value);
    }

}
