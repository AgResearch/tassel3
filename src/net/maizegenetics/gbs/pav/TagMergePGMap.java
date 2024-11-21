/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class TagMergePGMap extends TagGWASPGMap {
    boolean[] ifJMap;
	byte[] jChr;
	int[] jPos;
	double[] jBinomP;
	int[] jSigSNPNumBestChr;
	byte[] familyNum;
	long[] familyCode;
    
    
    public TagMergePGMap (String tagGWASPGMapFileS, String tagJointPGMapFileS) {
        super.readTxtFile(tagGWASPGMapFileS);
        TagJointPGMap jpgmap = new TagJointPGMap(tagJointPGMapFileS, false);
        this.mergeJointPGMap(jpgmap);
    }

    public TagMergePGMap (String tagMergePGMapFileS) {
        this.readTxtFile(tagMergePGMapFileS);
    }
    
    public void mkGWASQualityFile (double lCutoff, int outnum, String GWASQualityFileS) {
        try {
            int cnt = 0;
            BufferedWriter bw = new BufferedWriter (new FileWriter(GWASQualityFileS), 65536);
            bw.write("PChr\tPPos\tGChr\tGPos\tLRatio2");
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.pChr[i] < TagPMap.sChr || this.pChr[i] > TagPMap.eChr) continue;
                if (this.lRatio2[i] < lCutoff) continue;
                bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(gChr[i])+"\t"+String.valueOf(gPos[i])+"\t");
                bw.write(String.valueOf(this.lRatio2[i]));
                bw.newLine();
                cnt++;
                if (cnt == outnum) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
        System.out.println("GWAS quality file is written");
    }
    
    public void mkJointQualityFile (double pCutoff, int outnum, String JointQualityFileS) {
        try {
            int cnt = 0;
            BufferedWriter bw = new BufferedWriter (new FileWriter(JointQualityFileS), 65536);
            bw.write("PChr\tPPos\tJChr\tJPos\tJBinomP");
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.pChr[i] < TagPMap.sChr || this.pChr[i] > TagPMap.eChr) continue;
                if (this.jBinomP[i] > pCutoff) continue;
                bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                bw.write(String.valueOf(this.jBinomP[i]));
                bw.newLine();
                cnt++;
                if (cnt == outnum) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
        System.out.println("Joint linkage quality file is written");
    }
    
    public void mkResolutionFile (double lCutoff, double pCutoff, int outnum, String resolutionFileS) {
        try {
            int cnt = 0;
            BufferedWriter bw = new BufferedWriter (new FileWriter(resolutionFileS), 65536);
            bw.write("GPosLogDist\tJPosLogDist\tLRatio2\tJBinomP");
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.pChr[i] < TagPMap.sChr || this.pChr[i] > TagPMap.eChr) continue;
                if (this.lRatio2[i] < lCutoff) continue;
                if (this.jBinomP[i] > pCutoff) continue;
                if (!(pChr[i] == gChr[i] && pChr[i] == jChr[i])) continue;
                bw.write(String.valueOf(Math.abs(pPos[i]-gPos[i]))+"\t"+String.valueOf(Math.abs(pPos[i]-jPos[i]))+"\t"+String.valueOf(this.lRatio2[i])+"\t"+String.valueOf(this.jBinomP[i]));
                bw.newLine();
                cnt++;
                if (cnt == outnum) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
        System.out.println("Resolution file is written");
    }
    
    @Override
    public void readTxtFile (String tagMergePGMapFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(tagMergePGMapFileS), 65536);
            String temp;
            int tagNum = -1;
            while ((temp = br.readLine()) != null) tagNum++;
            this.iniMatrix(tagNum);
            br = new BufferedReader (new FileReader(tagMergePGMapFileS), 65536);
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
                ifJMap[i] = Boolean.valueOf(tem[19]);
                if (ifJMap[i]) {
                    jChr[i] = Byte.valueOf(tem[20]); jPos[i] = Integer.valueOf(tem[21]); jBinomP[i] = Double.valueOf(tem[22]);
                    jSigSNPNumBestChr[i] = Integer.valueOf(tem[23]); familyNum[i] = Byte.valueOf(tem[24]); familyCode[i] = Long.valueOf(tem[25]);
                }
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is read from " + tagMergePGMapFileS);
    }
    
    @Override
    public void iniMatrix (int tagNum) {
        super.iniMatrix(tagNum);
        ifJMap = new boolean[tagNum];
        jChr = new byte[tagNum];
        jPos = new int[tagNum];
        jBinomP = new double[tagNum];
        jSigSNPNumBestChr = new int[tagNum];
        familyNum = new byte[tagNum];
        familyCode = new long[tagNum];
    }
    
    public void mergeJointPGMap (TagJointPGMap jpgmap) {
        ifJMap = new boolean[this.getTagCount()];
        jChr = new byte[this.getTagCount()];
        jPos = new int[this.getTagCount()];
        jBinomP = new double[this.getTagCount()];
        jSigSNPNumBestChr = new int[this.getTagCount()];
        familyNum = new byte[this.getTagCount()];
        familyCode = new long[this.getTagCount()];
        for (int i = 0; i < this.getTagCount(); i++) {
            long[] t = this.getTag(i);
            int index = jpgmap.getTagIndex(t);
            if (index < 0) {
                ifJMap[i] = false;
            }
            else {
                ifJMap[i] = true;
                jChr[i] = jpgmap.jChr[index][0];
                jPos[i] = jpgmap.jPos[index][0];
                jBinomP[i] = jpgmap.jBinomP[index][0];
                jSigSNPNumBestChr[i] = jpgmap.jSigSNPNumChr[index][0];
                familyNum[i] = jpgmap.familyNum[index][0];
                familyCode[i] = jpgmap.familyCode[index][0];
            }
            
        }
    }
    
    @Override
    public void writeTxtFile (String tagGWASPGMapFileS) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "IfJMap\tJChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagGWASPGMapFileS), 65536);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\t");
                bw.write(String.valueOf(ifJMap[i])+"\t");
                if (ifJMap[i]) {
                    bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                    bw.write(String.valueOf(jBinomP[i]) +"\t"+String.valueOf(jSigSNPNumBestChr[i])+"\t"+String.valueOf(familyNum[i])+"\t");
                    bw.write(String.valueOf(familyCode[i]));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagGWASPGMapFileS);
    }
    
    public void writeTxtFile (String tagGWASPGMapFileS, boolean[] ifOut) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "IfJMap\tJChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagGWASPGMapFileS), 65536);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!(ifOut[i])) continue;
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\t");
                bw.write(String.valueOf(ifJMap[i])+"\t");
                if (ifJMap[i]) {
                    bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                    bw.write(String.valueOf(jBinomP[i]) +"\t"+String.valueOf(jSigSNPNumBestChr[i])+"\t"+String.valueOf(familyNum[i])+"\t");
                    bw.write(String.valueOf(familyCode[i]));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagGWASPGMapFileS);
    }
    
     public void writeTxtFile (String tagGWASPGMapFileS, boolean[] ifOut, int startIndex, int size) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "IfJMap\tJChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode";
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(tagGWASPGMapFileS), 65536);
            bw.write(header);
            bw.newLine();
            int cnt = 0;
            for (int i = startIndex; i < this.getTagCount(); i++) {
                if (!(ifOut[i])) continue;
                long[] t = new long[2];
                for (int j = 0; j < this.tagLengthInLong; j++) t[j] = tags[j][i];
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+String.valueOf(tagLength[i])+"\t");
                bw.write(String.valueOf(pChr[i]) +"\t"+String.valueOf(pPos[i])+"\t");
                bw.write(String.valueOf(ifMatch[i]) +"\t"+String.valueOf(ifPerfectMatch[i])+"\t"+String.valueOf(ifUniqueMatch[i])+"\t");
                bw.write(String.valueOf(gChr[i]) +"\t"+String.valueOf(gSite[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(gBinomP[i])+"\t");
                bw.write(String.valueOf(gSigSNPNum[i]) +"\t"+String.valueOf(tagTaxaCount[i])+"\t"+String.valueOf(sigChrNum[i])+"\t");
                bw.write(String.valueOf(lRatio2[i]) +"\t"+String.valueOf(lRatioM[i])+"\t"+String.valueOf(gSigSNPNumBestChr[i])+"\t");
                bw.write(String.valueOf(gSigSNPBCStartPos[i]) +"\t"+String.valueOf(gSigSNPBCEndPos[i])+"\t");
                bw.write(String.valueOf(ifJMap[i])+"\t");
                if (ifJMap[i]) {
                    bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t");
                    bw.write(String.valueOf(jBinomP[i]) +"\t"+String.valueOf(jSigSNPNumBestChr[i])+"\t"+String.valueOf(familyNum[i])+"\t");
                    bw.write(String.valueOf(familyCode[i]));
                }
                bw.newLine();
                cnt++;
                if (cnt == size) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is written to " + tagGWASPGMapFileS);
    }
}
