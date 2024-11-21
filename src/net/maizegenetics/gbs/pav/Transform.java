/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.util.Arrays;

/**
 *
 * @author Fei Lu
 */
public class Transform {
    long[] region;
    double[] rate;
    
    public Transform (String recombinationFileS) {
        Table t = new Table (recombinationFileS);
        region = new long[t.getRowNumber()];
        rate = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            region[i] = this.transformPos(Byte.valueOf(t.content[i][0]), Integer.valueOf(t.content[i][1]));
            rate[i] = Double.valueOf(t.content[i][2]);
        }
    }
    
    private long transformPos (byte chr, int pos) {
        long tPos = (long)chr * 1000000000 + pos;
        return tPos;
    }
    
    double getTagCount (int cnt) {
        return this.boxcox(cnt, 0.0202);
    }
    
    double getTagTaxaCount (int cnt) {
        return this.boxcox(cnt, 0.1);
    }
    
    double getGRecom (long tranPos) {
        int hit = Arrays.binarySearch(region, tranPos);
        if (hit < 0) hit = -hit - 2;
        return this.getGRecom(rate[hit]);
    }
    
    double getGRecom (double v) {
        if (v < 8.73E-06) v = 1E-07;
        return this.boxcox(v, 0.303);
    }
    
    double getGBinomP (double v) {
        if (v < 2.31E-308) v = 2.31E-308;
        return this.boxcox(v, 0.005);
    }
    
    double getLRatio2 (double v) {
        if (v > 300) v = 300;
        if (new Double(v).isNaN()) v = 300;
        return this.boxcox(v, 0.2);
    }
    
    double getLRatioM (double v) {
        if (v > 300) v = 300;
        if (new Double(v).isNaN()) v = 300; 
        return this.boxcox(v, 0);
    }
    
    double getGSigSNPNum (int cnt) {
        if (cnt == 0) cnt = 1;
        return this.boxcox(cnt, 0.22);
    }
    
    double getGSigSNPNumBC (int cnt) {
        if (cnt == 0) cnt = 1;
        return this.boxcox(cnt, -0.1);
    }
    
    double getGWidth (long cnt) {
        if (cnt == 0) cnt = 1;
        return this.boxcox(cnt, 0.1);
    }
    
    double getJRecom (long tranPos) {
        int hit = Arrays.binarySearch(region, tranPos);
        if (hit < 0) hit = -hit - 2;
        return this.getJRecom(rate[hit]);
    }
    
    double getJRecom (double v) {
        if (v < 8.73E-06) v = 1E-07;
        return this.boxcox(v, 0.3);
    }
    
    double getJBinomP (double v) {
        if (v < 2.31E-308) v = 2.31E-308;
        return this.boxcox(v, 0.15);
    }
    
    double getJSigSNPNumBC (int cnt) {
        return this.boxcox(cnt, 0.4);
    }
    
    double getFamilyNum (int cnt) {
        return this.boxcox(cnt, 0);
    }
    
    double getGDist (long cnt) {
        if (cnt == 0) cnt = 1;
        return Math.log10(cnt);
    }
    
    double getJDist (long cnt) {
        if (cnt == 0) cnt = 1;
        return Math.log10(cnt);
    }
    
    double getGJDist (long cnt) {
        if (cnt == 0) cnt = 1;
        return this.boxcox(cnt, 0.1);
    }
    
    private double boxcox (double y, double lambda) {
        if (lambda != 0) {
            return (Math.pow(y, lambda)-1)/lambda;
        }
        else {
            return Math.log(y);
        }
    }
}
