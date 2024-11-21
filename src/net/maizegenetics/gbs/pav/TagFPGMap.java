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
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Used to store tags with high resolution, which are filtered by model prediction (GJ,G and J)
 * Record the tags mapped by either GWAS or Joint, the fields not in the other mapping result are given -1.
 * @author Fei Lu
 */
public class TagFPGMap extends TagMPGMap {
    float[] predictionValue;
    
    
    public TagFPGMap (String tagMPGMapFileS, String tagGWASPGMapFileS, String tagJointPGMapFileS, String gjpgMapResPredictionValueFileS, String gpgMapResPredictionValueFileS, String jpgMapResPredictionValueFileS) {
        TagMPGMap mpgmap = new TagMPGMap (tagMPGMapFileS);
        TagGWASPGMap gpgmap = new TagGWASPGMap (tagGWASPGMapFileS);
        TagJointPGMap jpgmap = new TagJointPGMap (tagJointPGMapFileS, false);
        this.mergeHighResPGMap(mpgmap, gpgmap, jpgmap, gjpgMapResPredictionValueFileS, gpgMapResPredictionValueFileS,jpgMapResPredictionValueFileS);
    }
    
    
    public TagFPGMap (String tagFPGMapFileS) {
        this.readTxtFile(tagFPGMapFileS);
    }
    
    @Override
    public void iniMatrix (int tagNum) {
        super.iniMatrix(tagNum);
        predictionValue = new float[tagNum];
    }
    
    public void mkResolutionFile (int outnum, String resolutionFileS, boolean ifOnlyGWAS) {
        try {
            int cnt = 0;
            BufferedWriter bw = new BufferedWriter (new FileWriter(resolutionFileS), 65536);
            bw.write("LogPosDist");
            bw.newLine();
            int stop = 0;
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!this.getIfRefUniqueTag(i)) continue;
                long tgPos;
                if (this.getIfMapByGWAS(i)) {
                    tgPos = this.getTransformPos(gChr[i], gPos[i]);
                }
                else {
                    if (ifOnlyGWAS) continue;
                    tgPos = this.getTransformPos(jChr[i], jPos[i]);
                }
                long tpPos = this.getTransformPos(pChr[i], pPos[i]);
                long gd = Math.abs(tpPos-tgPos);
                if (gd < 65) continue;
                bw.write(String.valueOf(Math.log10(gd)));
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
    
    public void mkFPGMappingQuality (String outfileS, int outNum) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("pChr\tpPos\tgenChr\tgenPos");
            bw.newLine();
            int cnt = 0, stop = 0;
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.getIfRefUniqueTag(i)) {
                    if (this.getIfMapByGWAS(i)) {
                        bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(gChr[i])+"\t"+String.valueOf(gPos[i]));
                    }
                    else {
                        bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i]));
                    }
                    
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
            System.out.println("FPGMap mapping check is written to " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);;
        }
    }
    
    public void writePAVPosFileS (String PAVPosFileS, String tagCountFileS) {
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(PAVPosFileS), 65536);
            bw.write("Chr\tPos\tcount");
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = this.getTag(i);
                int hit = tc.getTagIndex(t);
                int cnt = tc.getReadCount(hit);
                if (this.getIfMapByGWAS(i)) {
                    bw.write(String.valueOf(gChr[i])+"\t"+String.valueOf(gPos[i])+"\t"+String.valueOf(cnt));
                }
                else {
                    bw.write(String.valueOf(jChr[i])+"\t"+String.valueOf(jPos[i])+"\t"+String.valueOf(cnt));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("PAVPos file is written at " + PAVPosFileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void writeAlignmentCheck2 (String samFileS, int distance, String alignFrequencyFileS) {
        System.out.println("Start reading sam file");
        int[] rightCount = new int[12];
        int[] alignCount = new int[12];
        for (int i = 0; i < rightCount.length; i++) rightCount[i] = 0;
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            ArrayList<String> tagAlnList = new ArrayList();
            ArrayList<Integer> tagIDList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                if (tagAlnList.isEmpty()) {
                    tagAlnList.add(temp);
                    String[] tem = temp.split("\t");
                    tagIDList.add(Integer.valueOf(tem[0]));
                }
                else {
                    String[] tem = temp.split("\t");
                    int tempID = Integer.valueOf(tem[0]);
                    if (tempID == tagIDList.get(0)) {
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                    else {
                        int index = getAlignmentIndex(tagAlnList, tagIDList.get(0), distance);
                        if (index == 0) {
                            alignCount[index]++;
                            rightCount[index]++;
                        }
                        else if (index == 11) {
                            alignCount[11]++;
                            for (int j = 0; j < tagIDList.size(); j++) {
                                alignCount[j+1]++;
                            }
                            rightCount[index]++;
                        }
                        else {
                            for (int j = 0; j < tagIDList.size(); j++) {
                                alignCount[j+1]++;
                            }
                            rightCount[index]++;
                        }
                        tagAlnList.clear();
                        tagIDList.clear();
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                }
            }
            int index = getAlignmentIndex(tagAlnList, tagIDList.get(0), distance);
            if (index == 0) {
                alignCount[index]++;
                rightCount[index]++;
            }
            else if (index == 11) {
                alignCount[11]++;
                for (int j = 0; j < tagIDList.size(); j++) {
                    alignCount[j+1]++;
                }
                rightCount[index]++;
            }
            else {
                for (int j = 0; j < tagIDList.size(); j++) {
                    alignCount[j+1]++;
                }
                rightCount[index]++;
            }
        }
        catch (Exception e) {System.out.println(e.toString());}
        System.out.println("Sam file read. Physical position of tags imported");
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(alignFrequencyFileS), 65536);
            int sum = 0;
            for (int i = 0; i < rightCount.length; i++) {
                bw.write(String.valueOf(i)+"Aln\t");
                sum += rightCount[i];
            }
            bw.write(String.valueOf(distance)+"bp");
            bw.newLine();
            for (int i = 0; i < rightCount.length; i++) {
                bw.write(String.valueOf((double)rightCount[i])+"\t");
            }
            bw.newLine();
            for (int i = 0; i < rightCount.length; i++) {
                bw.write(String.valueOf((double)alignCount[i])+"\t");
            }
            bw.newLine();
            for (int i = 0; i < rightCount.length; i++) {
                bw.write(String.valueOf(((double)rightCount[i]/alignCount[i]))+"\t");
            }
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println("AlignCheck file is written at " + alignFrequencyFileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void writeAlignmentCheck (String samFileS, int distance, String alignFrequencyFileS) {
        System.out.println("Start reading sam file");
        int[] count = new int[11];
        for (int i = 0; i < count.length; i++) count[i] = 0;
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            ArrayList<String> tagAlnList = new ArrayList();
            ArrayList<Integer> tagIDList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                if (tagAlnList.isEmpty()) {
                    tagAlnList.add(temp);
                    String[] tem = temp.split("\t");
                    tagIDList.add(Integer.valueOf(tem[0]));
                }
                else {
                    String[] tem = temp.split("\t");
                    int tempID = Integer.valueOf(tem[0]);
                    if (tempID == tagIDList.get(0)) {
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                    else {
                        int index = getAlignmentIndex(tagAlnList, tagIDList.get(0), distance);
                        count[index]++;
                        tagAlnList.clear();
                        tagIDList.clear();
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                }
            }
            int index = getAlignmentIndex(tagAlnList, tagIDList.get(0), distance);
            count[index]++;
        }
        catch (Exception e) {System.out.println(e.toString());}
        System.out.println("Sam file read. Physical position of tags imported");
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(alignFrequencyFileS), 65536);
            int sum = 0;
            for (int i = 0; i < count.length; i++) {
                bw.write(String.valueOf(i)+"Aln\t");
                sum += count[i];
            }
            bw.write(String.valueOf(distance)+"bp");
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                bw.write(String.valueOf((double)count[i]/sum)+"\t");
            }
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println("AlignCheck file is written at " + alignFrequencyFileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    private int getAlignmentIndex (ArrayList<String> tagAlnList, int tagId, int distance) {
        int index = tagId - 1;
        String[] temp = tagAlnList.get(0).split("\t");
        if (temp[2].startsWith("*")) {
            return 0;
        }
        for (int i = 0; i < tagAlnList.size(); i++) {
            temp = tagAlnList.get(i).split("\t");
            int chr = Integer.valueOf(temp[2]);
            int posi = Integer.valueOf(temp[3]);
            
            String c = temp[5]+"?";
            String[] tem = c.split("M");
            int match = 0;
            for (int j = 0; j < tem.length-1; j++) {
                tem[j] = tem[j].replaceFirst(".+\\D", "");
                match += Integer.valueOf(tem[j]);
            }
            if (this.getIfMapByGWAS(index)) {
                if (chr == this.gChr[index] && Math.abs(posi - this.gPos[index]) < distance) {
                    return i + 1;
                }
            }
            else {
                if (chr == this.jChr[index] && Math.abs(posi - this.jPos[index]) < distance) {
                    return i + 1;
                }
            }
        }
        return 11;
    }
    
    public void screenPrintProportionOfMappedTags (String samFileS) {
        boolean[] ifPAV = new boolean[this.getTagCount()];
        int[] count = new int[9];
        for (int i = 0; i < count.length; i++) count[i] = 0;
        System.out.println("Start reading sam file");
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            ArrayList<String> tagAlnList = new ArrayList();
            ArrayList<Integer> tagIDList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                if (tagAlnList.isEmpty()) {
                    tagAlnList.add(temp);
                    String[] tem = temp.split("\t");
                    tagIDList.add(Integer.valueOf(tem[0]));
                }
                else {
                    String[] tem = temp.split("\t");
                    int tempID = Integer.valueOf(tem[0]);
                    if (tempID == tagIDList.get(0)) {
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                    else {
                        //checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
                        count[this.getPortionCode(tagAlnList, tagIDList.get(0), ifPAV)]++;
                        tagAlnList.clear();
                        tagIDList.clear();
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                }
            }
            //checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
            count[this.getPortionCode(tagAlnList, tagIDList.get(0), ifPAV)]++;
        }
        catch (Exception e) {System.out.println(e.toString());}
        for (int i = 0; i < count.length; i++) System.out.println(count[i]);
        System.out.println("Sam file read. Physical position of tags imported");
    }
    
    public boolean[] getIfPAV (String samFileS) {
        boolean[] ifPAV = new boolean[this.getTagCount()];
        System.out.println("Start reading sam file");
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            ArrayList<String> tagAlnList = new ArrayList();
            ArrayList<Integer> tagIDList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                if (tagAlnList.isEmpty()) {
                    tagAlnList.add(temp);
                    String[] tem = temp.split("\t");
                    tagIDList.add(Integer.valueOf(tem[0]));
                }
                else {
                    String[] tem = temp.split("\t");
                    int tempID = Integer.valueOf(tem[0]);
                    if (tempID == tagIDList.get(0)) {
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                    else {
                        //checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
                        checkIfAbsentPAV(tagAlnList, tagIDList.get(0), ifPAV);
                        tagAlnList.clear();
                        tagIDList.clear();
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                }
            }
            //checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
            checkIfAbsentPAV(tagAlnList, tagIDList.get(0), ifPAV);
        }
        catch (Exception e) {System.out.println(e.toString());}
        System.out.println("Sam file read. Physical position of tags imported");
        return ifPAV;
    }
    
    private int getPortionCode (ArrayList<String> tagAlnList, int tagId, boolean[] ifPAV) {
        int distance = 10000000;
        int index = tagId - 1;
        String[] temp = tagAlnList.get(0).split("\t");
        if (temp[2].startsWith("*")) {
            return 0;
        }
        int min = 5;
        if (tagAlnList.size() < min) min = tagAlnList.size();
        int perfectB73Count = 0;
        for (int i = 0; i < min; i++) {
            temp = tagAlnList.get(i).split("\t");
            String testS = String.valueOf(this.getTagLength(index))+"M";
            if (temp[5].equals(testS) && tagAlnList.get(i).contains("NM:i:0")) {
                perfectB73Count++;
            }
        }
        
        for (int i = 0; i < min; i++) {
            temp = tagAlnList.get(i).split("\t");
            int chr = Integer.valueOf(temp[2]);
            int posi = Integer.valueOf(temp[3]);
            if (this.getIfMapByGWAS(index)) {
                if (chr == this.gChr[index] && Math.abs(posi - this.gPos[index]) < distance) {
                    if (perfectB73Count == 0) {
                        return 3; //NonB73Match tags agree with genetic position
                    }
                    else if (perfectB73Count == 1) {
                        return 1; //Unique B73 tags agree with genetic position
                    }
                    else {
                        return 2; //Multiple hits B73 tags agree with genetic position
                    }
                }
            }
            else {
                if (chr == this.jChr[index] && Math.abs(posi - this.jPos[index]) < distance) {
                    if (perfectB73Count == 0) {
                        return 3; //NonB73Match tags agree with genetic position
                    }
                    else if (perfectB73Count == 1) {
                        return 1; //Unique B73 tags agree with genetic position
                    }
                    else {
                        return 2; //Multiple hits B73 tags agree with genetic position
                    }
                }
            }
        }
        if (perfectB73Count == 0) {
            return 6; //NonB73Match tags don't agree with genetic position
        }
        else if (perfectB73Count == 1) {
            return 4; //Unique B73 tags don't agree with genetic position
        }
        else {
            return 5; //Multiple hits B73 tags don't agree with genetic position
        }
    }
    
    private void checkIfAbsentPAV (ArrayList<String> tagAlnList, int tagId, boolean[] ifPAV) {
        int distance = 10000000;
        int index = tagId - 1;
        String[] temp = tagAlnList.get(0).split("\t");
        if (temp[2].startsWith("*")) {
            ifPAV[index] = true;
            return;
        }
        ifPAV[index] = true;
        int min = 5;
        if (tagAlnList.size() < min) min = tagAlnList.size();
        for (int i = 0; i < min; i++) {
            temp = tagAlnList.get(i).split("\t");
            int chr = Integer.valueOf(temp[2]);
            int posi = Integer.valueOf(temp[3]);
            
            String c = temp[5]+"?";
            String[] tem = c.split("M");
            int match = 0;
            for (int j = 0; j < tem.length-1; j++) {
                tem[j] = tem[j].replaceFirst(".+\\D", "");
                match += Integer.valueOf(tem[j]);
            }
            if (this.getIfMapByGWAS(index)) {
                if (chr == this.gChr[index] && Math.abs(posi - this.gPos[index]) < distance) {
                    ifPAV[index] = false;
                    return;
                }
            }
            else {
                if (chr == this.jChr[index] && Math.abs(posi - this.jPos[index]) < distance) {
                    ifPAV[index] = false;
                    return;
                }
            }
        }
    }
    
    private void checkIfPAV (ArrayList<String> tagAlnList, int tagId, boolean[] ifPAV) {
        int distance = 500000;
        int minMatch = 40;
        int index = tagId - 1;
        String[] temp = tagAlnList.get(0).split("\t");
        if (temp[2].startsWith("*")) {
            ifPAV[index] = true;
            return;
        }
        ifPAV[index] = true;
        for (int i = 0; i < tagAlnList.size(); i++) {
            temp = tagAlnList.get(0).split("\t");
            int chr = Integer.valueOf(temp[2]);
            int posi = Integer.valueOf(temp[3]);
            
            String c = temp[5]+"?";
            String[] tem = c.split("M");
            int match = 0;
            for (int j = 0; j < tem.length-1; j++) {
                tem[j] = tem[j].replaceFirst(".+\\D", "");
                match += Integer.valueOf(tem[j]);
            }
            if (this.getIfMapByGWAS(index)) {
                if (chr == this.gChr[index] && Math.abs(posi - this.gPos[index]) < distance && match > minMatch) {
                    ifPAV[index] = false;
                    return;
                }
            }
            else {
                if (chr == this.jChr[index] && Math.abs(posi - this.jPos[index]) < distance && match > minMatch) {
                    ifPAV[index] = false;
                    return;
                }
            }
        }
    }
    
    public void writeFastaFileWithTagCountID (String tagCountFileS, String outfileS) {
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] temp = new long[this.tagLengthInLong];
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                int index = tc.getTagIndex(temp);
                String name = ">"+String.valueOf(index+1);
                bw.write(name);
                bw.newLine();
                bw.write(BaseEncoder.getSequenceFromLong(temp).substring(0, this.getTagLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Fasta file of tagPGMap output in " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
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
                predictionValue[i] = Float.valueOf(tem[25]);
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        System.out.println("TagGWASPGMap is read from " + tagMPGMapFileS);
    }
    
    @Override
    public void writeTxtFile (String tagMPGMapFileS) {
        this.writeTxtFile(tagMPGMapFileS, 0, this.getTagCount());
    }
    
    @Override
    public void writeTxtFile (String tagMPGMapFileS, int startIndex, int size) {
        String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tGChr\tSite\tGPos\tBinomP\tSigSNPNum\tTagTaxaCount\tSigChrNum\tLRatio2\tLRatioM\tSigSNPNumBC\tMinSigBCPos\tMaxSigBCPos\t"
                + "JChr\tJPos\tJBinomP\tJSigSNPNumBC\tFamilyNum\tFamilyCode\tpredictionValue";
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
                bw.write(String.valueOf(familyCode[i])+"\t"+String.valueOf(predictionValue[i]));
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
    
    public void writeFPGPaperFile (String outfileS, boolean[] ifPAV) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write(String.valueOf(this.tagLengthInLong));
            bw.newLine();
            bw.write(String.valueOf(this.getTagCount()));
            bw.newLine();
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tLog10(PredictedDistance)");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
                String pavS;
                if (ifPAV[i]) pavS = "1";
                else pavS = "0";
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                if (this.gChr[i] == -1) {
                    bw.write(String.valueOf(this.jChr[i])+"\t"+String.valueOf(this.jPos[i]));
                }
                else {
                    bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i]));
                }
                bw.write("\t"+pavS+"\t"+String.valueOf(predictionValue[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeTextTOGMFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            //bw.write(String.valueOf(this.tagLengthInLong));
            //bw.newLine();
            //bw.write(String.valueOf(this.getTagCount()));
            //bw.newLine();
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tPredictedDistance");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                if (this.gChr[i] == -1) {
                    bw.write(String.valueOf(this.jChr[i])+"\t"+String.valueOf(this.jPos[i]));
                }
                else {
                    bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i]));
                }
                bw.write("\t"+String.valueOf(Byte.MIN_VALUE)+"\t"+String.valueOf(this.predictionValue[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeFastaFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                String name = ">"+String.valueOf(i+1);
                bw.write(name);
                bw.newLine();
                long[] temp = new long[this.tagLengthInLong];
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp).substring(0, this.getTagLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Fasta file of tagPGMap output in " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void checkMappingQuality () {
        int B73Cnt = 0, cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same500k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            
            cnt++;
            if (this.getIfRefUniqueTag(i)) {
                B73Cnt++;
                if (this.getIfMapByGWAS(i)) {
                    if (!(this.pChr[i] == this.gChr[i])) continue;
                    sameChrCnt++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 100000) same100k++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 200000) same200k++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 500000) same500k++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 1000000) same1m++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 2000000) same2m++;
                    if (Math.abs(this.pPos[i] - this.gPos[i]) < 10000000) same10m++;
                }
                else {
                    if (!(this.pChr[i] == this.jChr[i])) continue;
                    sameChrCnt++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 100000) same100k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 200000) same200k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 500000) same500k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 1000000) same1m++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 2000000) same2m++;
                    if (Math.abs(this.pPos[i] - this.jPos[i]) < 10000000) same10m++;
                }
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
    
    private boolean getIfMapByGWAS (int index) {
        if (this.gChr[index] == -1) return false;
        return true;
    }
    
    
    private void mergeHighResPGMap (TagMPGMap mpgmap, TagGWASPGMap gpgmap,  TagJointPGMap jpgmap, String gjpgMapResPredictionValueFileS, String gpgMapResPredictionValueFileS, String jpgMapResPredictionValueFileS) {
        boolean[] ifInM = new boolean[mpgmap.getTagCount()];
        for (int i = 0; i < ifInM.length; i++) ifInM[i] = true;
        boolean[] ifInG = new boolean[gpgmap.getTagCount()];
        for (int i = 0; i < ifInG.length; i++) ifInG[i] = true;
        boolean[] ifInJ = new boolean[jpgmap.getTagCount()];
        for (int i = 0; i < ifInJ.length; i++) ifInJ[i] = true;
        for (int i = 0; i < ifInM.length; i++) {
            long[] tm = mpgmap.getTag(i);
            int hit = gpgmap.getTagIndex(tm);
            if (hit >= 0) ifInG[hit] = false;
            hit = jpgmap.getTagIndex(tm);
            if (hit >= 0) ifInJ[hit] = false;
        }
        int tagNum = this.getMergeSize(ifInM, ifInG, ifInJ);
        this.iniMatrix(tagNum);
        int index = this.mergeMPGMap(mpgmap, ifInM, 0, gjpgMapResPredictionValueFileS);
        index = this.mergeGWASPGMap(gpgmap, ifInG, index, gpgMapResPredictionValueFileS);
        this.mergeJointPGMap(jpgmap, ifInJ, index, jpgMapResPredictionValueFileS);
        this.sortByTag();
    }
    
    private int mergeJointPGMap (TagJointPGMap jpgmap, boolean[] ifInJ, int startIndex, String jpgMapResPredictionValueFileS) {
        Table t = new Table(jpgMapResPredictionValueFileS);
        int index = startIndex;
        for (int i = 0; i < jpgmap.getTagCount(); i++) {
            if (!ifInJ[i]) continue;
            this.importJointPGRecord(jpgmap, i, index, Float.valueOf(t.content[i][1]));
            index++;
        }
        return index;
    }
    
    private void importJointPGRecord (TagJointPGMap jpgmap, int i, int index, float predict) {
        long[] t = jpgmap.getTag(i);
        for (int j = 0; j < t.length; j++) tags[j][index] = t[j];
        tagLength[index] = (byte)jpgmap.getTagLength(i);
        pChr[index] = jpgmap.pChr[i];
        pPos[index] = jpgmap.pPos[i];
        ifMatch[index] = jpgmap.ifMatch[i];
        ifPerfectMatch[index] = jpgmap.ifPerfectMatch[i];
        ifUniqueMatch[index] = jpgmap.ifUniqueMatch[i];
        gChr[index] = -1;
        gSite[index] = -1;
        gPos[index] = -1;
        gBinomP[index] = -1;
        gSigSNPNum[index] = -1;
        tagTaxaCount[index] = -1;
        sigChrNum[index] = -1;
        lRatio2[index] = -1;
        lRatioM[index] = -1;
        gSigSNPNumBestChr[index] = -1;
        gSigSNPBCStartPos[index] = -1;
        gSigSNPBCEndPos[index] = -1;
        jChr[index] = jpgmap.jChr[i][0];
        jPos[index] = jpgmap.jPos[i][0];
        jBinomP[index] = jpgmap.jBinomP[i][0];
        jSigSNPNumBestChr[index] = jpgmap.jSigSNPNumChr[i][0];
        familyNum[index] = jpgmap.familyNum[i][0];
        familyCode[index] = jpgmap.familyCode[i][0];
        predictionValue[index] = predict;
    }
    
    private int mergeGWASPGMap (TagGWASPGMap gpgmap, boolean[] ifInG, int startIndex, String gpgMapResPredictionValueFileS) {
        Table t = new Table(gpgMapResPredictionValueFileS);
        int index = startIndex;
        for (int i = 0; i < gpgmap.getTagCount(); i++) {
            if (!ifInG[i]) continue;
            this.importJointPGRecord(gpgmap, i, index, Float.valueOf(t.content[i][1]));
            index++;
        }
        return index;
    }
    
    private void importJointPGRecord (TagGWASPGMap gpgmap, int i, int index, float predict) {
        long[] t = gpgmap.getTag(i);
        for (int j = 0; j < t.length; j++) tags[j][index] = t[j];
        tagLength[index] = (byte)gpgmap.getTagLength(i);
        pChr[index] = gpgmap.pChr[i];
        pPos[index] = gpgmap.pPos[i];
        ifMatch[index] = gpgmap.ifMatch[i];
        ifPerfectMatch[index] = gpgmap.ifPerfectMatch[i];
        ifUniqueMatch[index] = gpgmap.ifUniqueMatch[i];
        gChr[index] = gpgmap.gChr[i];
        gSite[index] = gpgmap.gSite[i];
        gPos[index] = gpgmap.gPos[i];
        gBinomP[index] = gpgmap.gBinomP[i];
        gSigSNPNum[index] = gpgmap.gSigSNPNum[i];
        tagTaxaCount[index] = gpgmap.tagTaxaCount[i];
        sigChrNum[index] = gpgmap.sigChrNum[i];
        lRatio2[index] = gpgmap.lRatio2[i];
        lRatioM[index] = gpgmap.lRatioM[i];
        gSigSNPNumBestChr[index] = gpgmap.gSigSNPNumBestChr[i];
        gSigSNPBCStartPos[index] = gpgmap.gSigSNPBCStartPos[i];
        gSigSNPBCEndPos[index] = gpgmap.gSigSNPBCEndPos[i];
        jChr[index] = -1;
        jPos[index] = -1;
        jBinomP[index] = -1;
        jSigSNPNumBestChr[index] = -1;
        familyNum[index] = -1;
        familyCode[index] = -1;
        predictionValue[index] = predict;
    }
    
    private int mergeMPGMap (TagMPGMap mpgmap, boolean[] ifInM, int startIndex, String gjpgMapResPredictionValueFileS) {
        Table t = new Table(gjpgMapResPredictionValueFileS);
        int index = startIndex;
        for (int i = 0; i < mpgmap.getTagCount(); i++) {
            if (!ifInM[i]) continue;
            this.importMPGRecord(mpgmap, i, index, Float.valueOf(t.content[i][1]));
            index++;
        }
        return index;
    }
    
    private void importMPGRecord (TagMPGMap mpgmap, int i, int index, float predict) {
        long[] t = mpgmap.getTag(i);
        for (int j = 0; j < t.length; j++) tags[j][index] = t[j];
        tagLength[index] = (byte)mpgmap.getTagLength(i);
        pChr[index] = mpgmap.pChr[i];
        pPos[index] = mpgmap.pPos[i];
        ifMatch[index] = mpgmap.ifMatch[i];
        ifPerfectMatch[index] = mpgmap.ifPerfectMatch[i];
        ifUniqueMatch[index] = mpgmap.ifUniqueMatch[i];
        gChr[index] = mpgmap.gChr[i];
        gSite[index] = mpgmap.gSite[i];
        gPos[index] = mpgmap.gPos[i];
        gBinomP[index] = mpgmap.gBinomP[i];
        gSigSNPNum[index] = mpgmap.gSigSNPNum[i];
        tagTaxaCount[index] = mpgmap.tagTaxaCount[i];
        sigChrNum[index] = mpgmap.sigChrNum[i];
        lRatio2[index] = mpgmap.lRatio2[i];
        lRatioM[index] = mpgmap.lRatioM[i];
        gSigSNPNumBestChr[index] = mpgmap.gSigSNPNumBestChr[i];
        gSigSNPBCStartPos[index] = mpgmap.gSigSNPBCStartPos[i];
        gSigSNPBCEndPos[index] = mpgmap.gSigSNPBCEndPos[i];
        jChr[index] = mpgmap.jChr[i];
        jPos[index] = mpgmap.jPos[i];
        jBinomP[index] = mpgmap.jBinomP[i];
        jSigSNPNumBestChr[index] = mpgmap.jSigSNPNumBestChr[i];
        familyNum[index] = mpgmap.familyNum[i];
        familyCode[index] = mpgmap.familyCode[i];
        predictionValue[index] = predict;
    }
    
    private int getMergeSize (boolean[] ifInM, boolean[] ifInG) {
        int cnt = 0;
        for (int i = 0; i < ifInM.length; i++) {
            if (ifInM[i]) cnt++;
        }
        for (int i = 0; i < ifInG.length; i++) {
            if (ifInG[i]) cnt++;
        }
        return cnt;
    }
    
    private int getMergeSize (boolean[] ifInM, boolean[] ifInG, boolean[] ifInJ) {
        int cnt = 0;
        for (int i = 0; i < ifInM.length; i++) {
            if (ifInM[i]) cnt++;
        }
        for (int i = 0; i < ifInG.length; i++) {
            if (ifInG[i]) cnt++;
        }
        for (int i = 0; i < ifInJ.length; i++) {
            if (ifInJ[i]) cnt++;
        }
        return cnt;
    }
    
    @Override
    public void swap(int index1, int index2) {
        super.swap(index1, index2);
        float temFloat;
        temFloat = predictionValue[index1]; predictionValue[index1] = predictionValue[index2]; predictionValue[index2] = temFloat;
    }
}
