package com.cloudinside.bio.VcfToolbox.test;

public class LogTest {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(LogTest.class);

    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        log.debug("hello");
        System.out.println("x");
        System.err.println("y");
    }

}
