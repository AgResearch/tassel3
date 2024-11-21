/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Only record the tags mapped by both Joint and GWAS.
 * @author Fei Lu
 */
public class TagMPGMap extends TagGWASPGMap {
    byte[] jChr;
	int[] jPos;
	double[] jBinomP;
	int[] jSigSNPNumBestChr;
	byte[] familyNum;
	long[] familyCode;
    
    public TagMPGMap () {}
    
    public TagMPGMap (String tagGWASPGMapFileS, String tagJointPGMapFileS) {
        TagGWASPGMap gpgmap = new TagGWASPGMap (tagGWASPGMapFileS);
        TagJointPGMap jpgmap = new TagJointPGMap(tagJointPGMapFileS, false);
        this.mergePGMap(gpgmap, jpgmap);
        this.sortByTag();
    }
    
    public TagMPGMap (String tagMPGMapFileS) {
        this.readTxtFile(tagMPGMapFileS);
    }
     
    @Override
    public void mkWekaFile (String tagCountFileS, String recombinationFileS, String wekaFileS) {
        TagCounts tc = new TagCounts(tagCountFileS, TagsByTaxa.FilePacking.Byte);
        Transform tr = new Transform(recombinationFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(wekaFileS), 65536);
            bw.write("@relation GJ\n\n");
            bw.write("@attribute TagCount numeric\n");
            bw.write("@attribute TagTaxaCount numeric\n");
            bw.write("@attribute GRecom numeric\n");
            bw.write("@attribute GBinomP numeric\n");
            bw.write("@attribute lRatio2 numeric\n");
            bw.write("@attribute lRatioM numeric\n");
            bw.write("@attribute GSigSNPNum numeric\n");
            bw.write("@attribute GSigSNPNumBC numeric\n");
            bw.write("@attribute GWidth numeric\n");
            bw.write("@attribute JRecom numeric\n");
            bw.write("@attribute JBinomP numeric\n");
            bw.write("@attribute JSigSNPNumBC numeric\n");
            bw.write("@attribute FamilyNum numeric\n");
            bw.write("@attribute GJDist numeric\n");
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
                long tgPos = this.transformPos(gChr[i], gPos[i]);
                long tjPos = this.transformPos(jChr[i], jPos[i]);
                long gDist;
                if (pChr[i] < 1) {
                    gDist = Long.MAX_VALUE;
                }
                else {
                    long tpPos = this.transformPos(pChr[i], pPos[i]);
                    gDist = Math.abs(tpPos-tgPos);
                } 
                long gjDist = Math.abs(tjPos-tgPos);
                bw.write(String.valueOf(tr.getGWidth(width))+",");
                bw.write(String.valueOf(tr.getJRecom(tjPos))+",");
                bw.write(String.valueOf(tr.getJBinomP(jBinomP[i]))+",");
                bw.write(String.valueOf(tr.getJSigSNPNumBC(jSigSNPNumBestChr[i]))+",");
                bw.write(String.valueOf(tr.getFamilyNum(familyNum[i]))+",");
                bw.write(String.valueOf(tr.getGJDist(gjDist)+","));
                bw.write(String.valueOf(tr.getGDist(gDist)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Weka file is created at " + wekaFileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void mkResolutionFile (double lCutoff, double pCutoff, int outnum, String resolutionFileS) {
        try {
            int cnt = 0;
            BufferedWriter bw = new BufferedWriter (new FileWriter(resolutionFileS), 65536);
            bw.write("LogGPosDist\tLogJPosDist\tLRatio2\tJBinomP");
            bw.newLine();
            int stop = 0;
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.lRatio2[i] < lCutoff) continue;
                if (this.jBinomP[i] > pCutoff) continue;
                if (!this.getIfRefUniqueTag(i)) continue;
                long tpPos = this.getTransformPos(pChr[i], pPos[i]);
                long tgPos = this.getTransformPos(gChr[i], gPos[i]);
                long tjPos = this.getTransformPos(jChr[i], jPos[i]);
                long gd = Math.abs(tpPos-tgPos);
                long jd = Math.abs(tpPos-tjPos);
                if (gd < 65 || jd < 65) continue;
                bw.write(String.valueOf(Math.log10(gd))+"\t"+String.valueOf(Math.log10(jd))+"\t"+String.valueOf(this.lRatio2[i])+"\t"+String.valueOf(this.jBinomP[i]));
                bw.newLine();
                cnt++;
                if (cnt == outnum) {
                    stop = i;
                    break;
                }
            }
            bw.flush();
            bw.close();
            System.out.println(String.valueOf((double)stop/this.getTagCount()) + " tags are scanned");
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
        System.out.println("Resolution file is written");
    }
    
    @Override
    public void writeTxtFile (String tagMPGMapFileS, boolean[] ifOut) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "JChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagMPGMapFileS), 65536);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!ifOut[i]) continue;
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\t");
                bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                bw.write(String.valueOf(jBinomP[i]) +"\t"+String.valueOf(jSigSNPNumBestChr[i])+"\t"+String.valueOf(familyNum[i])+"\t");
                bw.write(String.valueOf(familyCode[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagMPGMapFileS);
    }
    
    public void writeTxtFile (String tagMPGMapFileS, int startIndex, int size) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "JChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagMPGMapFileS), 65536);
            bw.write(header);
            bw.newLine();
            for (int i = startIndex; i < startIndex+size; i++) {
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\t");
                bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                bw.write(String.valueOf(jBinomP[i]) +"\t"+String.valueOf(jSigSNPNumBestChr[i])+"\t"+String.valueOf(familyNum[i])+"\t");
                bw.write(String.valueOf(familyCode[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagMPGMapFileS);
    }
    
    @Override
    public void writeTxtFile (String tagMPGMapFileS) {
        this.writeTxtFile(tagMPGMapFileS, 0, this.getTagCount());
    }
    
    @Override
    public void readTxtFile (String tagMPGMapFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(tagMPGMapFileS), 65536);
            String temp;
            int tagNum = -1;
            while ((temp = br.readLine()) != null) tagNum++;
            this.iniMatrix(tagNum);
            br = new BufferedReader (new FileReader(tagMPGMapFileS), 65536);
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
                jChr[i] = Byte.valueOf(tem[19]); jPos[i] = Integer.valueOf(tem[20]); jBinomP[i] = Double.valueOf(tem[21]);
                jSigSNPNumBestChr[i] = Integer.valueOf(tem[22]); familyNum[i] = Byte.valueOf(tem[23]); familyCode[i] = Long.valueOf(tem[24]);
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is read from " + tagMPGMapFileS);
    }
    
    private void mergePGMap (TagGWASPGMap gpgmap, TagJointPGMap jpgmap) {
        boolean[] ifShared = this.getIfShared(gpgmap, jpgmap);
        int tagNum = 0;
        for (int i = 0; i < ifShared.length; i++) {
            if (ifShared[i]) tagNum++;
        }
        this.iniMatrix(tagNum);
        int cnt = 0;
        for (int i = 0; i < ifShared.length; i++) {
            if (!ifShared[i]) continue;
            long[] t = gpgmap.getTag(i);
            for (int j = 0; j < this.getTagSizeInLong(); j++) {
                tags[j][cnt] = t[j];
            }
            tagLength[cnt] = (byte)gpgmap.getTagLength(i);
            pChr[cnt] = gpgmap.pChr[i];
            gSite[cnt] = gpgmap.gSite[i];
            pPos[cnt] = gpgmap.pPos[i];
            ifMatch[cnt] = gpgmap.ifMatch[i];
            ifPerfectMatch[cnt] = gpgmap.ifPerfectMatch[i];
            ifUniqueMatch[cnt] = gpgmap.ifUniqueMatch[i];
            gChr[cnt] = gpgmap.gChr[i];
            gPos[cnt] = gpgmap.gPos[i];
            gBinomP[cnt] = gpgmap.gBinomP[i];
            gSigSNPNum[cnt] = gpgmap.gSigSNPNum[i];
            tagTaxaCount[cnt] = gpgmap.tagTaxaCount[i];
            sigChrNum[cnt] = gpgmap.sigChrNum[i];
            lRatio2[cnt] = gpgmap.lRatio2[i];
            lRatioM[cnt] = gpgmap.lRatioM[i];
            gSigSNPNumBestChr[cnt] = gpgmap.gSigSNPNumBestChr[i];
            gSigSNPBCStartPos[cnt] = gpgmap.gSigSNPBCStartPos[i];
            gSigSNPBCEndPos[cnt] = gpgmap.gSigSNPBCEndPos[i];
            int hit = jpgmap.getTagIndex(t);
            jChr[cnt] = jpgmap.jChr[hit][0];
            jPos[cnt] = jpgmap.jPos[hit][0];
            jBinomP[cnt] = jpgmap.jBinomP[hit][0];
            jSigSNPNumBestChr[cnt] = jpgmap.jSigSNPNumChr[hit][0];
            familyNum[cnt] = jpgmap.familyNum[hit][0];
            familyCode[cnt] = jpgmap.familyCode[hit][0];
            cnt++;
        }
    }
    
    @Override
    public void iniMatrix (int tagNum) {
        super.iniMatrix(tagNum);
        jChr = new byte[tagNum];
        jPos = new int[tagNum];
        jBinomP = new double[tagNum];
        jSigSNPNumBestChr = new int[tagNum];
        familyNum = new byte[tagNum];
        familyCode = new long[tagNum];
    }
    
    private boolean[] getIfShared (TagGWASPGMap gpgmap, TagJointPGMap jpgmap) {
        boolean[] ifShare = new boolean[gpgmap.getTagCount()];
        for (int i = 0; i < gpgmap.getTagCount(); i++) {
            long[] t = gpgmap.getTag(i);
            int hit = jpgmap.getTagIndex(t);
            if (hit < 0) ifShare[i] = false;
            else ifShare[i] = true;
        }
        return ifShare;
    }
    
    @Override
    public void swap(int index1, int index2) {
        super.swap(index1, index2);
        byte temByte;
        temByte = jChr[index1]; jChr[index1] = jChr[index2]; jChr[index2] = temByte;
        int temInt;
        temInt = jPos[index1]; jPos[index1] = jPos[index2]; jPos[index2] = temInt;
        double temDouble;
        temDouble = jBinomP[index1]; jBinomP[index1] = jBinomP[index2]; jBinomP[index2] = temDouble;
        temInt = jSigSNPNumBestChr[index1]; jSigSNPNumBestChr[index1] = jSigSNPNumBestChr[index2]; jSigSNPNumBestChr[index2] = temInt;
        temByte = familyNum[index1]; familyNum[index1] = familyNum[index2]; familyNum[index2] = temByte;
        long tempLong;
        tempLong = familyCode[index1]; familyCode[index1] = familyCode[index2]; familyCode[index2] = tempLong;
    } 
}
