/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import cern.colt.GenericSorting;
import cern.jet.stat.Probability;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Tags don't align have the -1 at pChr and pPos. Chromosomes from 0 to 12. 
 * Tags align to multiple position were assigned the pChr and pPos of the best position
 * The field of gSite should be deleted, it's not useful
 * @author Fei Lu
 */
public class TagGWASPGMap extends AbstractTags {
    byte[] pChr; //physical chromosome
	int[] pPos; //physical pos by bowtie2
	boolean[] ifMatch;
	boolean[] ifPerfectMatch;
	boolean[] ifUniqueMatch;
    byte[] gChr; //genetic chromosome
	int[] gSite; //marker ID of this genetic position
	int[] gPos; //genetic position
    double[] gBinomP; //p-value of binomial test
    int[] gSigSNPNum; //number of significant SNP
	int[] tagTaxaCount; //how many taxa has this tag
	byte[] sigChrNum; //chromosome numbers with significant SNP
    double[] lRatio2; //likelihood ratio of the best chr verses the second best chr
    double[] lRatioM; //likelihood ratio of the best chr verses the median best chr
	int[] gSigSNPNumBestChr;
	int[] gSigSNPBCStartPos; //significant SNP start position on the best chromosome
    int[] gSigSNPBCEndPos;
    
    public TagGWASPGMap () {}
    
    public TagGWASPGMap (String tagPMapFileS, String tagGWASGMapFileS, String tagCountFileS) {
        TagPMap pMap = new TagPMap (tagPMapFileS);
		TagCounts tc = new TagCounts (tagCountFileS, TagsByTaxa.FilePacking.Byte);
        this.readGWASResult(tagGWASGMapFileS, tc, pMap);
        this.sortByTag();
    }
    
    public TagGWASPGMap (String tagGWASPGMapFileS) {
        readTxtFile(tagGWASPGMapFileS);
        this.sortByTag();
    }
    
    public boolean[] getIfHighResolution (String predictionFileS, double disCut) {
        boolean[] ifOut = new boolean[this.getTagCount()];
        try {
            BufferedReader br = new BufferedReader (new FileReader(predictionFileS), 65536);
            br.readLine();
            String temp;
            for (int i = 0; i < ifOut.length; i++) {
                temp = br.readLine();
                String[] tem = temp.split("\t");
                if (Double.valueOf(tem[2]) <= disCut) ifOut[i] = true;
                else ifOut[i] = false;
            }
        }
        catch(Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        return ifOut;
    }
    
    public void mkWekaFile (String tagCountFileS, String recombinationFileS, String predictionFileS) {
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        Transform tr = new Transform(recombinationFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(predictionFileS), 65536);
            bw.write("@relation G\n\n");
            bw.write("@attribute TagCount numeric\n");
            bw.write("@attribute TagTaxaCount numeric\n");
            bw.write("@attribute GRecom numeric\n");
            bw.write("@attribute GBinomP numeric\n");
            bw.write("@attribute lRatio2 numeric\n"); 
            bw.write("@attribute lRatioM numeric\n");
            bw.write("@attribute GSigSNPNum numeric\n");
            bw.write("@attribute GSigSNPNumBC numeric\n");
            bw.write("@attribute GWidth numeric\n");
            bw.write("@attribute GDist numeric\n\n");
            bw.write("@data\n");
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = this.getTag(i);
                int hit = tc.getTagIndex(t);
                bw.write(String.valueOf(tr.getTagCount(tc.getReadCount(hit)))+",");
                bw.write(String.valueOf(tr.getTagTaxaCount(this.tagTaxaCount[i]))+",");
                bw.write(String.valueOf(tr.getGRecom(this.getTransformPos(gChr[i], gPos[i])))+",");
                bw.write(String.valueOf(tr.getGBinomP(this.gBinomP[i]))+",");
                bw.write(String.valueOf(tr.getLRatio2(this.lRatio2[i]))+",");
                bw.write(String.valueOf(tr.getLRatioM(this.lRatioM[i]))+",");
                bw.write(String.valueOf(tr.getGSigSNPNum(this.gSigSNPNum[i]))+",");
                bw.write(String.valueOf(tr.getGSigSNPNumBC(this.gSigSNPNumBestChr[i]))+",");
                long start = this.transformPos(gChr[i], gSigSNPBCStartPos[i]);
                long end = this.transformPos(gChr[i], gSigSNPBCEndPos[i]);
                long width = Math.abs(start-end)+1;
                long gDist;
                if (pChr[i] < 1) {
                    gDist = Long.MAX_VALUE;
                }
                else {
                    long tpPos = this.transformPos(pChr[i], pPos[i]);
                    long tgPos = this.transformPos(gChr[i], gPos[i]);
                    gDist = Math.abs(tpPos-tgPos);
                }
                bw.write(String.valueOf(tr.getGWidth(width))+",");
                bw.write(String.valueOf(tr.getGDist(gDist)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    long transformPos (byte chr, int pos) {
        long tPos = (long)chr * 1000000000 + pos;
        return tPos;
    }
     
    public void mkGMappingQuality (String outfileS, double pThresh, double lTresh, int outNum) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("pChr\tpPos\tgenChr\tgenPos");
            bw.newLine();
            int cnt = 0, stop = 0;
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.gBinomP[i] > pThresh) continue;
                if (this.lRatio2[i] < lTresh) continue;
                if (this.getIfRefUniqueTag(i)) {
                    bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(gChr[i])+"\t"+String.valueOf(gPos[i]));
                    bw.newLine();
                    cnt++;
                }
                if (cnt == outNum) {
                    stop = i;
                    break;
                }
            }
            bw.flush();
            bw.close();
            System.out.println(String.valueOf((double)stop/this.getTagCount())+" is scanned");
            System.out.println("GWAS mapping check is written to " + outfileS);
            
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void checkGMappingQuality (double pThresh, double lTresh) {
        int B73Cnt = 0, cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same500k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.gBinomP[i] > pThresh) continue;
            if (this.lRatio2[i] < lTresh) continue;
            cnt++;
            if (this.getIfRefUniqueTag(i)) {
                B73Cnt++;
                if (!(this.pChr[i] == this.gChr[i])) continue;
                sameChrCnt++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 100000) same100k++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 200000) same200k++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 500000) same500k++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 1000000) same1m++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 2000000) same2m++;
                if (Math.abs(this.pPos[i] - this.gPos[i]) < 10000000) same10m++;
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under the threshhold");
        System.out.println(B73Cnt+ " unique B73 tags are checked");
        System.out.println((double)sameChrCnt/B73Cnt + " (persentage) tags has agreement on genetic Chr and physical Chr");
        System.out.println((double)same10m/B73Cnt + " are in 10M region");
        System.out.println((double)same2m/B73Cnt + " are in 2M region");
        System.out.println((double)same1m/B73Cnt + " are in 1M region");
        System.out.println((double)same500k/B73Cnt + " are in 500k region");
        System.out.println((double)same200k/B73Cnt + " are in 200k region");
        System.out.println((double)same100k/B73Cnt + " are in 100k region");
    }
    
    public boolean getIfNonRefMatchTag (int index) {
        if (ifMatch[index] == true && this.ifPerfectMatch[index] == false && pChr[index] >= TagPMap.sChr && pChr[index] <= TagPMap.eChr) return true;
        return false;
    }
    
    public boolean getIfRefTag (int index) {
        if (ifMatch[index] == true && this.ifPerfectMatch[index] == true && pChr[index] >= TagPMap.sChr && pChr[index] <= TagPMap.eChr) return true;
        return false;
    }
    
    public boolean[] getIfRefTag() {
        boolean[] ifRef = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.getIfRefTag(i)) {
                ifRef[i] = true;
                cnt++;
            }
            else ifRef[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " uniquely aligned ref tags");
        return ifRef;
    }
    
    public boolean[] getFilterCross (boolean[] f1, boolean[] f2) {
        boolean[] fCross = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (f1[i] == true && f2[i] == true) {
                fCross[i] = true;
                cnt++;
            }
            else fCross[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " passed the cross filter");
        return fCross;
    }
    
    public boolean getIfRefUniqueTag (int index) {
        if (ifMatch[index] == true && ifPerfectMatch[index] == true && ifUniqueMatch[index] == true && pChr[index] >= TagPMap.sChr && pChr[index] <= TagPMap.eChr) return true;
        return false;
    }
    
    public boolean[] getIfRefUniqueTag () {
        boolean[] ifRef = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.getIfRefUniqueTag(i)) {
                ifRef[i] = true;
                cnt++;
            }
            else ifRef[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " uniquely aligned ref tags");
        return ifRef;
    }
    
    public boolean[] getIfPassLRatio2 (double lCutoff) {
        boolean[] ifRef = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.lRatio2[i] > lCutoff) {
                ifRef[i] = true;
                cnt++;
            }
            else ifRef[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " tags passed lRatio2 cutoff of " + String.valueOf(lCutoff));
        return ifRef;
    }
    
    long getTransformPos (byte chr, int pos) {
        long tPos = (long)chr * 1000000000 + pos;
        return tPos;
    }
    
    void mkCombinedP (TagJointPGMap jpgmap, String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
            bw.write("GWASP\tJointP\tCombinedX2\n");
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = this.getTag(i);
                int hit = jpgmap.getTagIndex(t);
                if (hit < 0) continue;
                double x2 = -2 * (Math.log(jpgmap.jBinomP[hit][0])+Math.log(this.gBinomP[i]));
                bw.write(String.valueOf(this.gBinomP[i])+"\t"+String.valueOf(jpgmap.jBinomP[hit][0])+"\t"+String.valueOf(1-Probability.chiSquare(4, x2)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
    }
    
    public void writeTxtFile (String tagGWASPGMapFileS) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\n";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagGWASPGMapFileS), 65536);
            bw.write(header);
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagGWASPGMapFileS);
    }
    
    public void writeTxtFile (String tagGWASPGMapFileS, boolean[] ifout) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\n";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagGWASPGMapFileS), 65536);
            bw.write(header);
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!ifout[i]) continue;
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagGWASPGMapFileS);
    }
    
    public void readTxtFile (String tagGWASPGMapFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(tagGWASPGMapFileS), 65536);
            String temp;
            int tagNum = -1;
            while ((temp = br.readLine()) != null) tagNum++;
            this.iniMatrix(tagNum);
            br = new BufferedReader (new FileReader(tagGWASPGMapFileS), 65536);
            br.readLine();
            for (int i = 0; i < tagNum; i++) {
                temp = br.readLine();
                String[] tem = temp.split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(tem[0]);
                for (int j = 0; j < t.length; j++) tags[j][i] = t[j];
                tagLength[i] = Byte.valueOf(tem[1]);
                pChr[i] = Byte.valueOf(tem[2]);
                pPos[i] = Integer.valueOf(tem[3]);
                ifMatch[i] = Boolean.valueOf(tem[4]); ifPerfectMatch[i] = Boolean.valueOf(tem[5]); ifUniqueMatch[i] = Boolean.valueOf(tem[6]);
                gChr[i] = Byte.valueOf(tem[7]); gSite[i] = Integer.valueOf(tem[8]); gPos[i] = Integer.valueOf(tem[9]);
                gBinomP[i] = Double.valueOf(tem[10]); gSigSNPNum[i] = Integer.valueOf(tem[11]); tagTaxaCount[i] = Integer.valueOf(tem[12]);
                sigChrNum[i] = Byte.valueOf(tem[13]); lRatio2[i] = Double.valueOf(tem[14]); lRatioM[i] = Double.valueOf(tem[15]);
                gSigSNPNumBestChr[i] = Integer.valueOf(tem[16]); gSigSNPBCStartPos[i] = Integer.valueOf(tem[17]); gSigSNPBCEndPos[i] = Integer.valueOf(tem[18]);
            }
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is read from " + tagGWASPGMapFileS);
    }
    
    private void readGWASResult (String tagGWASGMapFileS, TagCounts tc, TagPMap pMap) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(tagGWASGMapFileS), 65536);
            String temp;
            int tagNum = -1;
            while ((temp = br.readLine()) != null) tagNum++;
            this.iniMatrix(tagNum);
            br = new BufferedReader (new FileReader(tagGWASGMapFileS), 65536);
            temp = br.readLine();
            for (int i = 0; i < tagNum; i++) {
                temp = br.readLine();
                String[] tem = temp.split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(tem[0]);
                for (int j = 0; j < t.length; j++) tags[j][i] = t[j];
                int hit = tc.getTagIndex(t);
                tagLength[i] = (byte)tc.getTagLength(hit);
                hit = pMap.getTagIndex(t);
                if (hit < 0) {
                    pChr[i] = -1;
                    pPos[i] = -1;
                    ifMatch[i] = false;
                    ifPerfectMatch[i] = false;
                    ifUniqueMatch[i] = false;
                }
                else {
                    ifMatch[i] = true;
                    if (pMap.hitNum[hit] == 1) {
                        ifUniqueMatch[i] = true;
                    }
                    else {
                        ifUniqueMatch[i] = false;
                    }
                    ifPerfectMatch[i] = pMap.ifPerfectMatch[hit][0];
                    pChr[i] = pMap.pChr[hit][0];
                    pPos[i] = pMap.pPos[hit][0];
                }
                gChr[i] = Byte.valueOf(tem[5]);
                gSite[i] = Integer.valueOf(tem[6]);
                gPos[i] = Integer.valueOf(tem[7]);
                gBinomP[i] = Double.valueOf(tem[8]);
                gSigSNPNum[i] = Integer.valueOf(tem[9]);
                tagTaxaCount[i] = Integer.valueOf(tem[10]);
                sigChrNum[i] = Byte.valueOf(tem[11]);
                lRatio2[i] = Double.valueOf(tem[12]);
                lRatioM[i] = Double.valueOf(tem[13]);
                gSigSNPNumBestChr[i] = Integer.valueOf(tem[14]);
                gSigSNPBCStartPos[i] = Integer.valueOf(tem[15]);
                gSigSNPBCEndPos[i] = Integer.valueOf(tem[16]);
            }
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Tag GWAS result is read and physical postition is incorporated");
    }
   
    public void iniMatrix (int tagNum) {
        tagLengthInLong = 2;
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        pChr = new byte[tagNum];
        pPos = new int[tagNum];
        ifMatch = new boolean[tagNum];
        ifPerfectMatch = new boolean[tagNum];
        ifUniqueMatch = new boolean[tagNum];
        gChr = new byte[tagNum];
        gSite = new int[tagNum];
        gPos = new int[tagNum];
        gBinomP = new double[tagNum];
        gSigSNPNum = new int[tagNum];
        tagTaxaCount = new int[tagNum];
        sigChrNum = new byte[tagNum];
        lRatio2 = new double[tagNum];
        lRatioM = new double[tagNum];
        gSigSNPNumBestChr = new int[tagNum];
        gSigSNPBCStartPos = new int[tagNum];
        gSigSNPBCEndPos = new int[tagNum];
        System.out.println("Matrix of "+String.valueOf(tagNum)+" tags is initialized in TagGWASPGMap");
    }
    
    public void sortByTag() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
		System.out.println("TagPGMap file is sorted by tags");
    }
     
    @Override
    public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp=tags[i][index1];
            tags[i][index1]=tags[i][index2];
            tags[i][index2]=temp;
        }
        byte temByte;
        temByte = tagLength[index1]; tagLength[index1] = tagLength[index2]; tagLength[index2] = temByte;
        temByte = pChr[index1]; pChr[index1] = pChr[index2]; pChr[index2] = temByte;
		int temInt;
		temInt = pPos[index1]; pPos[index1] = pPos[index2]; pPos[index2] = temInt;
        boolean temBoo;
        temBoo = ifMatch[index1]; ifMatch[index1] = ifMatch[index2]; ifMatch[index2] = temBoo;
        temBoo = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = temBoo;
        temBoo = ifUniqueMatch[index1]; ifUniqueMatch[index1] = ifUniqueMatch[index2]; ifUniqueMatch[index2] = temBoo;
        temByte = gChr[index1]; gChr[index1] = gChr[index2]; gChr[index2] = temByte;
        temInt = gSite[index1]; gSite[index1] = gSite[index2]; gSite[index2] = temInt;
		temInt = gPos[index1]; gPos[index1] = gPos[index2]; gPos[index2] = temInt;
        double temDouble;
        temDouble = gBinomP[index1]; gBinomP[index1] = gBinomP[index2]; gBinomP[index2] = temDouble;
        temInt = gSigSNPNum[index1]; gSigSNPNum[index1] = gSigSNPNum[index2]; gSigSNPNum[index2] = temInt;
        temInt = tagTaxaCount[index1]; tagTaxaCount[index1] = tagTaxaCount[index2]; tagTaxaCount[index2] = temInt;
        temByte = sigChrNum[index1]; sigChrNum[index1] = sigChrNum[index2]; sigChrNum[index2] = temByte;
        temDouble = lRatio2[index1]; lRatio2[index1] = lRatio2[index2]; lRatio2[index2] = temDouble;
        temDouble = lRatioM[index1]; lRatioM[index1] = lRatioM[index2]; lRatioM[index2] = temDouble;
        temInt = gSigSNPNumBestChr[index1]; gSigSNPNumBestChr[index1] = gSigSNPNumBestChr[index2]; gSigSNPNumBestChr[index2] = temInt;
        temInt = gSigSNPBCStartPos[index1]; gSigSNPBCStartPos[index1] = gSigSNPBCStartPos[index2]; gSigSNPBCStartPos[index2] = temInt;
        temInt = gSigSNPBCEndPos[index1]; gSigSNPBCEndPos[index1] = gSigSNPBCEndPos[index2]; gSigSNPBCEndPos[index2] = temInt; 
    } 
}
