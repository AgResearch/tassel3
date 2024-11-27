/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import org.apache.commons.math.stat.ranking.NaturalRanking;

/**
 *
 * @author Fei Lu
 */
public class Train {
    TagMPGMap mpgmp;
    int distCut = 2000000;
    String header = "";
    long[] PPos;
    long[] GPos;
    long[] JPos;
    double[] tagCount;
    double[] tagTaxaCount;
    int[] hits;
    double[] GRecom;
    double[] GBinomP;
    double[] lRatio2;
    double[] lRatioM;
    double[] GSigSNPNum;
    byte[] GSigChrNum;
    double[] GSigSNPNumBC;
    double[] GWidth;
    double[] JRecom;
    double[] JBinomP;
    double[] JSigSNPNumBC;
    double[] familyNum;
    String[] totalClass;
    String[] GClass;
    String[] JClass;
    double[] GDist;
    double[] JDist;
    double[] GJDist;
    
    public Train (String tagMPGMapFileS) {
        this.addMergePGMap(tagMPGMapFileS);
    }
    
    public void writeTrainingFile (String trainingFileS) {
        String deli = ",";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(trainingFileS), 65536);
            String header = "PPos\tGPos\tJPos\tTagCount\tTagTaxaCount\tHits\t" +
                    "GRecom\tGBinomP\tlRatio2\tlRatioM\tGSigSNPNum\tGSigChrNum\tGSigSNPNumBC\tGWidth\t" +
                    "JRecom\tJBinomP\tJSigSNPNumBC\tFamilyNum\tTotalClass\tGClass\tJClass\tGDist\tJDist\tGJDist";
            header = header.replaceAll("\t", deli);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < this.getInstanceCount(); i++) {
                bw.write(String.valueOf(PPos[i])+deli+String.valueOf(GPos[i])+deli+String.valueOf(JPos[i])+deli);
                bw.write(String.valueOf(tagCount[i])+deli+String.valueOf(tagTaxaCount[i])+deli+String.valueOf(hits[i])+deli);
                bw.write(String.valueOf(GRecom[i])+deli+String.valueOf(GBinomP[i])+deli+String.valueOf(lRatio2[i])+deli);
                bw.write(String.valueOf(lRatioM[i])+deli+String.valueOf(GSigSNPNum[i])+deli+String.valueOf(GSigChrNum[i])+deli);
                bw.write(String.valueOf(GSigSNPNumBC[i])+deli+String.valueOf(GWidth[i])+deli+String.valueOf(JRecom[i])+deli);
                bw.write(String.valueOf(JBinomP[i])+deli+String.valueOf(JSigSNPNumBC[i])+deli+String.valueOf(familyNum[i])+deli+String.valueOf(totalClass[i])+deli);
                bw.write(String.valueOf(GClass[i])+deli+String.valueOf(JClass[i])+deli);
                bw.write(String.valueOf(GDist[i])+deli+String.valueOf(JDist[i])+deli+String.valueOf(GJDist[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
    }
    
    private double boxcox (double y, double lambda) {
        if (lambda != 0) {
            return (Math.pow(y, lambda)-1)/lambda;
        }
        else {
            return Math.log(y);
        }
    }
    
    public void boxcoxTransform () {
        for (int i = 0; i < this.getInstanceCount(); i++) {
            tagCount[i] = this.boxcox(tagCount[i], (double)0.0202);
            tagTaxaCount[i] = this.boxcox(tagTaxaCount[i], (double)0.1);
            GRecom[i] = this.boxcox(GRecom[i], (double)0.3030);
            GBinomP[i] = this.boxcox(GBinomP[i], (double)0.005);
            lRatio2[i] = this.boxcox(lRatio2[i], (double)0.2);
            lRatioM[i] = this.boxcox(lRatioM[i], (double)0.0);
            GSigSNPNum[i] = this.boxcox(GSigSNPNum[i], (double)0.22);
            GSigSNPNumBC[i] = this.boxcox(GSigSNPNumBC[i], (double)-0.1);
            GWidth[i] = this.boxcox(GWidth[i], (double)0.1);
            JRecom[i] = this.boxcox(JRecom[i], (double)0.3);
            JBinomP[i] = this.boxcox(JBinomP[i], (double)0.15);
            JSigSNPNumBC[i] = this.boxcox(JSigSNPNumBC[i], (double)0.4);
            familyNum[i] = this.boxcox(familyNum[i], 0);
            GDist[i] = Math.log10(GDist[i]);
            JDist[i] = Math.log10(JDist[i]);
            //JDist[i] = this.boxcox(JDist[i], (double)0.1); 
            GJDist[i] = this.boxcox(GJDist[i], (double)0.1);    
        }
    }
    
    public void rankTransform () {
        NaturalRanking ranking = new NaturalRanking();
        tagCount = ranking.rank(tagCount);
        tagTaxaCount = ranking.rank(tagTaxaCount);
        GRecom = ranking.rank(GRecom);
        GBinomP = ranking.rank(GBinomP);
        lRatio2 = ranking.rank(lRatio2);
        lRatioM = ranking.rank(lRatioM);
        GSigSNPNum = ranking.rank(GSigSNPNum);
        GSigSNPNumBC = ranking.rank(GSigSNPNumBC);
        GWidth = ranking.rank(GWidth);
        JRecom = ranking.rank(JRecom);
        JBinomP = ranking.rank(JBinomP);
        JSigSNPNumBC = ranking.rank(JSigSNPNumBC);
        for (int i = 0; i < this.getInstanceCount(); i++) {
            GDist[i] = Math.log10(GDist[i]);
            JDist[i] = Math.log10(JDist[i]);
        }
        GJDist = ranking.rank(GJDist); 
    }
    
    public void logTransform () {
        for (int i = 0; i < this.getInstanceCount(); i++) {
            tagCount[i] = Math.log10(tagCount[i]);
            tagTaxaCount[i] = Math.log10(tagTaxaCount[i]);
            if (GRecom[i] < 8.73E-06) GRecom[i] = 1E-07;
            GRecom[i] = -Math.log10(GRecom[i]);
            if (GBinomP[i] < 2.31E-308) GBinomP[i] = 2.31E-308;
            GBinomP[i] = -Math.log10(GBinomP[i]);
            if (lRatio2[i] > 300) lRatio2[i] = 300;
            if (new Double(lRatio2[i]).isNaN()) lRatio2[i] = 300;
            lRatio2[i] = Math.log10(lRatio2[i]);
            if (lRatioM[i] > 300) lRatioM[i] = 300;
            if (new Double(lRatioM[i]).isNaN()) lRatioM[i] = 300;
            lRatioM[i] = Math.log10(lRatioM[i]);
            if (GSigSNPNum[i] == 0) GSigSNPNum[i] = 1;
            GSigSNPNum[i] = Math.log10(GSigSNPNum[i]);
            if (GSigSNPNumBC[i] == 0) GSigSNPNumBC[i] = 1;
            GSigSNPNumBC[i] = Math.log10(GSigSNPNumBC[i]);
            GWidth[i] = Math.log10(GWidth[i]); //kind of skewed
            if (JRecom[i] < 8.73E-06) JRecom[i] = 1E-07;
            JRecom[i] = -Math.log10(JRecom[i]);
            if (JBinomP[i] < 2.31E-308) JBinomP[i] = 2.31E-308;
            JBinomP[i] = Math.log10(-Math.log10(JBinomP[i]));
            JSigSNPNumBC[i] = Math.log10(JSigSNPNumBC[i]); //kind of skewed
            if (GDist[i] == 0) GDist[i] = 1;
            if (JDist[i] == 0) JDist[i] = 1;
            if (GJDist[i] == 0) GJDist[i] = 1;
            GDist[i] = Math.log10(GDist[i]);
            JDist[i] = Math.log10(JDist[i]);
            GJDist[i] = Math.log10(GJDist[i]);
        }
    } 
    
    public void preTransform () {
         for (int i = 0; i < this.getInstanceCount(); i++) {
            if (GRecom[i] < 8.73E-06) GRecom[i] = 1E-07;
            if (GBinomP[i] < 2.31E-308) GBinomP[i] = 2.31E-308;
            if (lRatio2[i] > 300) lRatio2[i] = 300;
            if (new Double(lRatio2[i]).isNaN()) lRatio2[i] = 300;
            if (lRatioM[i] > 300) lRatioM[i] = 300;
            if (new Double(lRatioM[i]).isNaN()) lRatioM[i] = 300;
            if (GSigSNPNum[i] == 0) GSigSNPNum[i] = 1;
            if (GSigSNPNumBC[i] == 0) GSigSNPNumBC[i] = 1;
            if (GWidth[i] == 0) GWidth[i] = 1; //kind of skewed
            if (JRecom[i] < 8.73E-06) JRecom[i] = 1E-07;
            if (JBinomP[i] < 2.31E-308) JBinomP[i] = 2.31E-308;
            if (GDist[i] == 0) GDist[i] = 1;
            if (JDist[i] == 0) JDist[i] = 1;
            if (GJDist[i] == 0) GJDist[i] = 1;
        }
    }
    
    public void addRecombination (String recombinationFileS) {
        Table t = new Table (recombinationFileS);
        long[] region = new long[t.getRowNumber()];
        double[] rate = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            region[i] = this.transformPos(Byte.valueOf(t.content[i][0]), Integer.valueOf(t.content[i][1]));
            rate[i] = Double.valueOf(t.content[i][2]);
        }
        for (int i = 0; i < this.getInstanceCount(); i++) {
            int hit = Arrays.binarySearch(region, GPos[i]);
            if (hit < 0) hit = -hit - 2;
            GRecom[i] = rate[hit];
            hit = Arrays.binarySearch(region, JPos[i]);
            if (hit < 0) hit = -hit - 2;
            JRecom[i] = rate[hit];
        }
    }
    
    public void addTagCountHits (String tagCountFileS,String tagHitFileS) {
        TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Byte);
        int hit[] = new int[tc.getTagCount()];
        try { 
            BufferedReader br = new BufferedReader (new FileReader(tagHitFileS), 65536);
            int cnt = -1;
            String temp;
            while ((temp = br.readLine())!= null) cnt++;
            if (cnt != tc.getTagCount()) {
                System.err.println("tagHitFile doesn't have the right tag number");
                System.exit(1);
            }
            br = new BufferedReader (new FileReader(tagHitFileS), 65536);
            br.readLine();
            while ((temp = br.readLine())!= null) {
                String[] tem = temp.split("\t");
                int index = Integer.valueOf(tem[0]) -1;
                hit[index] = Integer.valueOf(tem[1]);
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
        for (int i = 0; i < mpgmp.getTagCount(); i++) {
            long[] t = mpgmp.getTag(i);
            int index = tc.getTagIndex(t);
            tagCount[i] = tc.getReadCount(index);
            hits[i] = hit[index];
        }
    }
    
    public void addMergePGMap (String tagMergePGMapFileS) {
        mpgmp = new TagMPGMap(tagMergePGMapFileS);
        this.iniMatrix();
        GBinomP = mpgmp.gBinomP;
        lRatio2 = mpgmp.lRatio2;
        lRatioM = mpgmp.lRatioM;
        GSigChrNum = mpgmp.sigChrNum;
        JBinomP = mpgmp.jBinomP;
        for (int i = 0; i < this.getInstanceCount(); i++) {
            tagTaxaCount[i] = mpgmp.tagTaxaCount[i];
            GSigSNPNum[i] = mpgmp.gSigSNPNum[i];
            GSigSNPNumBC[i] = mpgmp.gSigSNPNumBestChr[i];
            JSigSNPNumBC[i] = mpgmp.jSigSNPNumBestChr[i];
            familyNum[i] = mpgmp.familyNum[i];
        }
        
        for (int i = 0; i < this.getInstanceCount(); i++) {
            PPos[i] = this.transformPos(mpgmp.pChr[i], mpgmp.pPos[i]);
            GPos[i] = this.transformPos(mpgmp.gChr[i], mpgmp.gPos[i]);
            JPos[i] = this.transformPos(mpgmp.jChr[i], mpgmp.jPos[i]);
            GWidth[i] = Math.abs(mpgmp.gSigSNPBCStartPos[i] - mpgmp.gSigSNPBCEndPos[i]) + 1;
            long disG = Math.abs(PPos[i]-GPos[i]);
            long disJ = Math.abs(PPos[i]-JPos[i]);
            if (disG > distCut && disJ > distCut) {
                totalClass[i] = "W";
            }
            else {
                if (disG <= disJ) totalClass[i] = "G";
                else totalClass[i] = "J";
            }
            if (disG > distCut) GClass[i] = "W";
            else GClass[i] = "G";
            if (mpgmp.pChr[i] == mpgmp.jChr[i]) JClass[i] = "J";
            else JClass[i] = "W";
            GDist[i] = disG;
            JDist[i] = disJ;
            GJDist[i] = Math.abs(JPos[i]-GPos[i]);
        }
    }
    
    private long transformPos (byte chr, int pos) {
        long tPos = (long)chr * 1000000000 + pos;
        return tPos;
    }
    
    
    
    public void iniMatrix () {
        tagTaxaCount = new double[mpgmp.getTagCount()];
        GBinomP = mpgmp.gBinomP;
        lRatio2 = mpgmp.lRatio2;
        lRatioM = mpgmp.lRatioM;
        GSigSNPNum = new double[mpgmp.getTagCount()];
        GSigChrNum = mpgmp.sigChrNum;
        GSigSNPNumBC = new double[mpgmp.getTagCount()];
        JBinomP = mpgmp.jBinomP;
        JSigSNPNumBC = new double[mpgmp.getTagCount()];
        familyNum = new double[mpgmp.getTagCount()];
        PPos = new long[mpgmp.getTagCount()];
        GPos = new long[mpgmp.getTagCount()];
        JPos = new long[mpgmp.getTagCount()];
        tagCount = new double[mpgmp.getTagCount()];
        hits = new int[mpgmp.getTagCount()];
        GRecom = new double[mpgmp.getTagCount()];
        GWidth = new double[mpgmp.getTagCount()];
        JRecom = new double[mpgmp.getTagCount()];
        totalClass =new String[mpgmp.getTagCount()];
        GClass = new String[mpgmp.getTagCount()];
        JClass = new String[mpgmp.getTagCount()];
        GDist = new double[mpgmp.getTagCount()];
        JDist = new double[mpgmp.getTagCount()];
        GJDist = new double[mpgmp.getTagCount()];
    }
    
    public int getInstanceCount () {
        return PPos.length; 
    }
}
