/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.colt.GenericSorting;
import cern.colt.function.IntComparator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Tags don't align have the -1 at pChr and pPos. Chromosomes from 0 to 12. 
 * Tags align to multiple position were assigned the pChr and pPos of the best position
 * @author Fei Lu
 */
public class TagJointPGMap extends TagJMap {
	byte[] pChr;
	int[] pPos;
	boolean[] ifMatch;
	boolean[] ifPerfectMatch; // when there are multiple hits, it means if the first hit is perfect
	boolean[] ifUniqueMatch; // if it has only one hit

	public TagJointPGMap (String tagPMapFileS, String tagGMapFileS, String tagCountFileS) {
		super.readTxtFile(tagGMapFileS);
		TagPMap pMap = new TagPMap (tagPMapFileS);
		TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Byte);
		this.mergePMap(pMap, tc);
        this.orderFamilyGroupByP();
        this.sortByTag();
	}

	public TagJointPGMap (String tagPGMapFileS, boolean ifBinary) {
		if (ifBinary) this.readBinaryFile(tagPGMapFileS);
		else this.readTxtFile(tagPGMapFileS);
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
            System.out.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        return ifOut;
    }
    
    public void mkWekaFile (String tagCountFileS, String recombinationFileS, String wekaFileS) {
        TagCounts tc = new TagCounts(tagCountFileS, TagsByTaxa.FilePacking.Byte);
        Transform tr = new Transform(recombinationFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(wekaFileS), 65536);
            bw.write("@relation J\n\n");
            bw.write("@attribute TagCount numeric\n");
            bw.write("@attribute JRecom numeric\n");
            bw.write("@attribute JBinomP numeric\n");
            bw.write("@attribute JSigSNPNumBC numeric\n");
            bw.write("@attribute FamilyNum numeric\n");
            bw.write("@attribute JDist numeric\n\n");
            bw.write("@data\n");
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] t = this.getTag(i);
                int hit = tc.getTagIndex(t);
                bw.write(String.valueOf(tr.getTagCount(tc.getReadCount(hit)))+",");
         
                

                long tpPos = this.transformPos(pChr[i], pPos[i]);
                long tjPos = this.transformPos(jChr[i][0], jPos[i][0]);
                long jDist;
                if (pChr[i] < 1) {
                    jDist = Long.MAX_VALUE;
                }
                else {
                    jDist = Math.abs(tpPos-tjPos);
                }
                bw.write(String.valueOf(tr.getJRecom(tjPos))+",");
                bw.write(String.valueOf(tr.getJBinomP(jBinomP[i][0]))+",");
                bw.write(String.valueOf(tr.getJSigSNPNumBC(jSigSNPNumChr[i][0]))+",");
                bw.write(String.valueOf(tr.getFamilyNum(familyNum[i][0]))+",");
                bw.write(String.valueOf(tr.getJDist(jDist)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Weka file is created at " + wekaFileS);
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
    
	public void mergePMap (TagPMap pMap, TagCounts tc) {
		tagLength = new byte[this.getTagCount()];
		pChr = new byte[this.getTagCount()];
		pPos = new int[this.getTagCount()];
		ifMatch = new boolean[this.getTagCount()];
		ifPerfectMatch = new boolean[this.getTagCount()];
		ifUniqueMatch = new boolean[this.getTagCount()];
		for (int i = 0; i < this.getTagCount(); i++) {
			int index = pMap.getTagIndex(this.getTag(i));
			if (index < 0) {
				pChr[i] = -1;
				pPos[i] = -1;
				ifMatch[i] = false;
				ifPerfectMatch[i] = false;
				ifUniqueMatch[i] = false;
			}
			else {
				ifMatch[i] = true;
				if (pMap.hitNum[index] == 1) {
					ifUniqueMatch[i] = true;
				}
				else {
					ifUniqueMatch[i] = false;
				}
				ifPerfectMatch[i] = pMap.ifPerfectMatch[index][0];
				pChr[i] = pMap.pChr[index][0];
				pPos[i] = pMap.pPos[index][0];
			}
			int hit = tc.getTagIndex(this.getTag(i));
			tagLength[i] = (byte)tc.getTagLength(hit);
		}
		System.out.println("TagGMap and TagPMap is merged");
	}

    public void orderFamilyGroupByReference () {
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jNum[i] == 1) continue;
            int[] index = new int[this.jNum[i]];
            for (int j = 0; j < index.length; j++) {
                index[j] = j;
            }
            
            double[] p = new double[this.jNum[i]];
            
            
            for (int j = 0; j < p.length - 1; j++) {
                for (int k = j + 1; k < p.length; k++) {
                    if (p[j] > p[k]) {
                        double tempP = p[j];
                        int tempIndex = index[j];
                        p[j] = p[k]; p[k] = tempP;
                        index[j] = index[k]; index[k] = tempIndex;
                    }
                }
            }
            this.orderFamilyGroupByIndex(index, i);
        }
        System.out.println("Family groups are ordered by p value");
    }
    
    public void orderFamilyGroupByP () {
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jNum[i] == 1) continue;            
            if (this.pChr[i] == this.jChr[i][0]) continue;
            int[] index = new int[this.familyGroupNum];
            index[0] = 1; index[1] = 0;
            this.orderFamilyGroupByIndex(index, i);
        }
        System.out.println("Family groups are ordered by p value");
    }
    
    private void orderFamilyGroupByIndex (int[] familyIndex, int tagIndex) {
        byte[] tempGChr = new byte[familyIndex.length];
        int[] tempSite = new int[familyIndex.length];
        int[] tempGPos = new int[familyIndex.length];
        double[] tempP = new double[familyIndex.length];
        int[] tempChrSig = new int[familyIndex.length];
        double[] tempLikelihood = new double[familyIndex.length];
        byte[] tempFamilyNum = new byte[familyIndex.length];
        long[] tempFamilyCode = new long[familyIndex.length];
        for (int i = 0; i < familyIndex.length; i++) {
            tempGChr[i] = this.jChr[tagIndex][familyIndex[i]];
            tempSite[i] = this.site[tagIndex][familyIndex[i]];
            tempGPos[i] = this.jPos[tagIndex][familyIndex[i]];
            tempP[i] = this.jBinomP[tagIndex][familyIndex[i]];
            tempChrSig[i] = this.jSigSNPNumChr[tagIndex][familyIndex[i]];
            tempLikelihood[i] = this.likelyhood[tagIndex][familyIndex[i]];
            tempFamilyNum[i] = this.familyNum[tagIndex][familyIndex[i]];
            tempFamilyCode[i] = this.familyCode[tagIndex][familyIndex[i]];
        }
        this.jChr[tagIndex] = tempGChr; 
        this.site[tagIndex] = tempSite;
        this.jPos[tagIndex] = tempGPos;
        this.jBinomP[tagIndex] = tempP;
        this.jSigSNPNumChr[tagIndex] = tempChrSig;
        this.likelyhood[tagIndex] = tempLikelihood;
        this.familyNum[tagIndex] = tempFamilyNum;
        this.familyCode[tagIndex] = tempFamilyCode;
        
    }
    
    public boolean[] getReverseFilter (boolean[] f1) {
        boolean[] f2 = new boolean[f1.length];
        int cnt = 0;
        for (int i = 0; i < f2.length; i++) {
            if (f1[i] == true) f2[i] = false;
            else {
                f2[i] = true;
                cnt++;
            }
        }
        System.out.println(String.valueOf(cnt) + " are in the reversed filter");
        return f2;
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
    
    public boolean[] getIfUnderThresholdFirstFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jBinomP[i][0] > pThresh) ifUnderThresh[i] = false;
            else {
                ifUnderThresh[i] = true;
                cnt++;
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in the first family group");
        return ifUnderThresh;
    }
    
    public boolean[] getIfUnderThresholdEitherFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            ifUnderThresh[i] = false;
            for (int j = 0; j < this.jNum[i]; j++) {
                if (this.jBinomP[i][0] <= pThresh) {
                    ifUnderThresh[i] = true;
                    cnt++;
                    break;
                }
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in either family group");
        return ifUnderThresh;
    }
    
    public boolean[] getIfUnderThresholdAllFamilyGroup (double pThresh) {
        boolean[] ifUnderThresh = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            ifUnderThresh[i] = true;
            int c = 0;
            for (int j = 0; j < this.jNum[i]; j++) {
                if (this.jBinomP[i][0] <= pThresh) {
                    c++;
                    
                }
            }
            if (c == this.jNum[i]) {
                ifUnderThresh[i] = true;
                cnt++;
            }
            else ifUnderThresh[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " tags are under threshold in all family groups");
        return ifUnderThresh;
    }
    
    public boolean[] getIfSingleMapping (double pThresh) {
        boolean[] ifSingle = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jNum[i] == 1) {
                if (this.jBinomP[i][0] <= pThresh) {
                    ifSingle[i] = true;
                    cnt++;
                }
                else ifSingle[i] = false;
            }
            else {
                int c = 0;
                for (int j = 0; j < this.jNum[i]; j++) {
                    if (this.jBinomP[i][j] <= pThresh) c++;
                }
                if (c == 1) {
                    ifSingle[i] = true;
                    cnt++;
                }
                else ifSingle[i] = false;
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are uniquely mapped");
        return ifSingle;
    }
    
    public boolean[] getIfMultipleMapping (double pThresh) {
        boolean[] ifMultiple = new boolean[this.getTagCount()];
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jNum[i] > 1) {
                int c = 0;
                for (int j = 0; j < this.jNum[i]; j++) {
                   if (this.jBinomP[i][j] <= pThresh) c++;
                }
                if (c > 1) {
                    ifMultiple[i] = true;
                    cnt++;
                }
                else ifMultiple[i] = false;
            }
            else ifMultiple[i] = false;
        }
        System.out.println(String.valueOf(cnt) + " tags are multiplly mapped");
        return ifMultiple;
    }
    
    public boolean getIfRefTag (int index) {
        if (ifMatch[index] == true && this.ifPerfectMatch[index] == true) return true;
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
    
     public boolean getIfRefUniqueTag (int index) {
        if (ifMatch[index] == true && ifPerfectMatch[index] == true && ifUniqueMatch[index] == true) return true;
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
    
    public boolean[] getIfMatchUniqueNonGene (String samFileS) {
        boolean[] ifMatch = new boolean[this.getTagCount()];
        System.out.println("Start reading sam file");
        int[] indexCount = new int[this.getTagCount()];
        boolean[] ifGene = new boolean[this.getTagCount()];
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                String[] tem = temp.substring(0, 70).split("\t");
                int index = Integer.valueOf(tem[0])-1;
                indexCount[index]++;
                if (tem[2].equals("*")) ifGene[index] = true;
                else ifGene[index] = false;    
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < ifMatch.length; i++) {
            if (indexCount[i] > 1) ifMatch[i] = false;
            else {
                if (ifGene[i]) ifMatch[i] = true;
                else ifMatch[i] = false;
            }
        }
        return ifMatch;
    }
    
    public boolean[] getIfMatchUniqueGene (String samFileS) {
        boolean[] ifMatch = new boolean[this.getTagCount()];
        System.out.println("Start reading sam file");
        int[] indexCount = new int[this.getTagCount()];
        boolean[] ifGene = new boolean[this.getTagCount()];
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                String[] tem = temp.substring(0, 70).split("\t");
                int index = Integer.valueOf(tem[0])-1;
                indexCount[index]++;
                if (tem[2].equals("*")) ifGene[index] = false;
                else ifGene[index] = true;    
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < ifMatch.length; i++) {
            if (indexCount[i] > 1) ifMatch[i] = false;
            else {
                if (ifGene[i]) ifMatch[i] = true;
                else ifMatch[i] = false;
            }
        }
        return ifMatch;
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
						checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
						tagAlnList.clear();
						tagIDList.clear();
						tagAlnList.add(temp);
						tagIDList.add(Integer.valueOf(tem[0]));
					}
				}
			}
			checkIfPAV(tagAlnList, tagIDList.get(0), ifPAV);
		}
		catch (Exception e) {System.out.println(e.toString());}
		System.out.println("Sam file read. Physical position of tags imported");
		return ifPAV;
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
			
			if (chr == this.jChr[index][0] && Math.abs(posi - this.jPos[index][0]) < distance && match > minMatch) {
				ifPAV[index] = false;
				return;
			}
		}
	}
    
    public void checkJMappingQualityMultiGroup (double pThresh, int familyGroupIndex) {
        int cnt = 0, B73Cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.jNum[i] == 1) continue;
            boolean flag = false;
            for (int j = 0; j < this.jNum[i]; j++) {
                if (this.jBinomP[i][j] > pThresh) {
                    flag = true;
                    continue;
                }
            }
            if (flag) continue;
            cnt++;
            if (this.getIfRefUniqueTag(i)) {
                B73Cnt++;
                if (!(this.pChr[i] == this.jChr[i][familyGroupIndex])) continue;
                sameChrCnt++;
                if (Math.abs(this.pPos[i] - this.jPos[i][familyGroupIndex]) < 50000) same100k++;
                if (Math.abs(this.pPos[i] - this.jPos[i][familyGroupIndex]) < 100000) same200k++;
                if (Math.abs(this.pPos[i] - this.jPos[i][familyGroupIndex]) < 500000) same1m++;
                if (Math.abs(this.pPos[i] - this.jPos[i][familyGroupIndex]) < 1000000) same2m++;
                if (Math.abs(this.pPos[i] - this.jPos[i][familyGroupIndex]) < 5000000) same10m++;
            }
        }
        System.out.println(String.valueOf(cnt) + " multiple mapped tags are under the threshhold of " + String.valueOf(pThresh) + " at each family group");
        System.out.println(sameChrCnt + "/"+ B73Cnt+ " tags are checked");
        System.out.println((double)sameChrCnt/B73Cnt + " (persentage) tags has agreement on genetic Chr and physical Chr");
        System.out.println((double)same10m/B73Cnt + " are in 10M region");
        System.out.println((double)same2m/B73Cnt + " are in 2M region");
        System.out.println((double)same1m/B73Cnt + " are in 1M region");
        System.out.println((double)same200k/B73Cnt + " are in 200k region");
        System.out.println((double)same100k/B73Cnt + " are in 100k region");
    }
    
    public void mkJMappingQuality (String outfileS, double pThresh, int outNum) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("pChr\tpPos\tgenChr\tgenPos");
            bw.newLine();
            int cnt = 0, stop = 0;
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.jBinomP[i][0] > pThresh) continue;
                if (this.getIfRefUniqueTag(i)) {
                    bw.write(String.valueOf(pChr[i])+"\t"+String.valueOf(pPos[i])+"\t"+String.valueOf(jChr[i][0])+"\t"+String.valueOf(jPos[i][0]));
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
            System.out.println("Joint mapping check is written to " + outfileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);;
        }
    }
    
    public void checkJMappingQuality (double pThresh, boolean ifFirstFamilyGroup, boolean ifSingleMapTag) {
        int B73Cnt = 0, cnt = 0;
        int sameChrCnt = 0;
        int same100k = 0, same200k = 0, same500k = 0, same1m = 0, same2m = 0, same10m = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (ifSingleMapTag) {
                if (this.jNum[i] != 1) continue;
            }
            if (ifFirstFamilyGroup) {
                if (this.jBinomP[i][0] > pThresh) continue;
                cnt++;
                if (this.getIfRefUniqueTag(i)) {
                    B73Cnt++;
                    if (!(this.pChr[i] == this.jChr[i][0])) continue;
                    sameChrCnt++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 100000) same100k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 200000) same200k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 500000) same500k++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 1000000) same1m++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 2000000) same2m++;
                    if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 10000000) same10m++;
                }
            }
            else {
                boolean flag = true;
                for (int j = 0; j < this.jNum[i]; j++) {
                    if (this.jBinomP[i][j] > pThresh) continue;
                    if (flag) {
                        cnt++;
                        flag = false;
                    }
                    if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) B73Cnt++;
                    break;
                }
                for (int j = 0; j < this.jNum[i]; j++) {
                    if (this.jBinomP[i][j] > pThresh) continue;
                    if (this.ifMatch[i] == true && this.ifPerfectMatch[i] && this.ifUniqueMatch[i] == true) {
                        if (!(this.pChr[i] == this.jChr[i][j])) continue;
                        sameChrCnt++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 100000) same100k++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 200000) same200k++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 500000) same500k++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 1000000) same1m++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 2000000) same2m++;
                        if (Math.abs(this.pPos[i] - this.jPos[i][0]) < 10000000) same10m++;
                        break;
                    }
                }
            }
        }
        System.out.println(String.valueOf(cnt) + " tags are under the threshhold of " + String.valueOf(pThresh));
        System.out.println(sameChrCnt + "/"+ B73Cnt+ " tags are checked");
        System.out.println((double)sameChrCnt/B73Cnt + " (persentage) tags has agreement on genetic Chr and physical Chr");
        System.out.println((double)same10m/B73Cnt + " are in 10M region");
        System.out.println((double)same2m/B73Cnt + " are in 2M region");
        System.out.println((double)same1m/B73Cnt + " are in 1M region");
        System.out.println((double)same500k/B73Cnt + " are in 500k region");
        System.out.println((double)same200k/B73Cnt + " are in 200k region");
        System.out.println((double)same100k/B73Cnt + " are in 100k region");
    }

	@Override
	public void readBinaryFile (String infileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
			int tagNum = dis.readInt();
			this.tagLengthInLong = dis.readInt();
			this.iniMatrix(tagNum);
			for (int i = 0; i < tagNum; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					tags[j][i] = dis.readLong();
				}
				tagLength[i] = dis.readByte();
				pChr[i] = dis.readByte();
				pPos[i] = dis.readInt();
				ifMatch[i] = dis.readBoolean();
				ifPerfectMatch[i] = dis.readBoolean();
				ifUniqueMatch[i] = dis.readBoolean();
				tagTaxaCount[i] = dis.readInt();
				jNum[i] = dis.readByte();
				this.iniSubMatrix(jNum[i], i);
				for (int j = 0; j < jNum[i]; j++) {
					jChr[i][j] = dis.readByte();
					site[i][j] = dis.readInt();
					jPos[i][j] = dis.readInt();
					jBinomP[i][j] = dis.readDouble();
					jSigSNPNumChr[i][j] = dis.readInt();
					likelyhood[i][j] = dis.readDouble();
					familyNum[i][j] = dis.readByte();
					familyCode[i][j] = dis.readLong();
				}
			}
			System.out.println("Binary tagPGMap file is read from " + infileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
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
			System.out.println("Fasta file of tagPGMap output");
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	@Override
	public void writeBinaryFile (String outfileS) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
			dos.writeInt(this.getTagCount());
			dos.writeInt(this.tagLengthInLong);
			for (int i = 0; i < this.getTagCount(); i++) {
				for (int j = 0; j < this.tagLengthInLong; j++){
					dos.writeLong(tags[j][i]);
				}
				dos.writeByte(tagLength[i]);
				dos.writeByte(pChr[i]);
				dos.writeInt(pPos[i]);
				dos.writeBoolean(ifMatch[i]);
				dos.writeBoolean(ifPerfectMatch[i]);
				dos.writeBoolean(ifUniqueMatch[i]);
				dos.writeInt(tagTaxaCount[i]);
				dos.writeByte(jNum[i]);
				for (int j = 0; j < jNum[i]; j++) {
					dos.writeByte(jChr[i][j]);
					dos.writeInt(site[i][j]);
					dos.writeInt(jPos[i][j]);
					dos.writeDouble(jBinomP[i][j]);
					dos.writeInt(jSigSNPNumChr[i][j]);
					dos.writeDouble(likelyhood[i][j]);
					dos.writeByte(familyNum[i][j]);
					dos.writeLong(familyCode[i][j]);
				}
			}
			dos.flush();
			dos.close();
			System.out.println("Binary tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	private void iniMatrix (int tagNum) {
		tags = new long[tagLengthInLong][tagNum];
		tagLength = new byte[tagNum];
		pChr = new byte[tagNum];
		pPos = new int[tagNum];
		ifMatch = new boolean[tagNum];
		ifPerfectMatch = new boolean[tagNum];
		ifUniqueMatch = new boolean[tagNum];
		tagTaxaCount = new int[tagNum];
		jNum = new byte[tagNum];
		jChr = new byte[tagNum][];
		site = new int[tagNum][];
		jPos = new int[tagNum][];
		jBinomP = new double[tagNum][];
		jSigSNPNumChr = new int[tagNum][];
		likelyhood = new double[tagNum][];
		familyNum = new byte[tagNum][];
		familyCode = new long[tagNum][];
		System.out.println("Matrix of tagPGMap is initialized with " + tagNum + " tags");
	}

	private void iniSubMatrix (int groupNum, int index) {
		jChr[index] = new byte[groupNum];
		site[index] = new int[groupNum];
		jPos[index] = new int[groupNum];
		jBinomP[index] = new double[groupNum];
		jSigSNPNumChr[index] = new int[groupNum];
		likelyhood[index] = new double[groupNum];
		familyNum[index] = new byte[groupNum];
		familyCode[index] = new long[groupNum];
	}

	@Override
	public void readTxtFile (String infileS) {
		int tagNum = this.getSize(infileS);
		this.iniMatrix(tagNum);
		try {
			BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
			String temp = br.readLine();
			for (int i = 0; i < tagNum; i++) {
				temp = br.readLine();
				this.readTxtRecord(temp, i);
			}
			System.out.println("Txt tagPGMap file is read from " + infileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	private void readTxtRecord (String line, int index) {
		String[] temp = line.split("\\s+");
		long[] temTag = BaseEncoder.getLongArrayFromSeq(temp[0]);
		for (int i = 0; i < this.tagLengthInLong; i++) {
			tags[i][index] = temTag[i];
		}
		tagLength[index] = Byte.valueOf(temp[1]);
		pChr[index] = Byte.valueOf(temp[2]);
		pPos[index] = Integer.valueOf(temp[3]);
		ifMatch[index] = Boolean.valueOf(temp[4]);
		ifPerfectMatch[index] = Boolean.valueOf(temp[5]);
		ifUniqueMatch[index] = Boolean.valueOf(temp[6]);
		tagTaxaCount[index] = Integer.valueOf(temp[7]);
		jNum[index] = Byte.valueOf(temp[8]);
		this.iniSubMatrix(jNum[index], index);
		for (int i = 0; i < jNum[index]; i++) {
			int j = i*8+9;
			jChr[index][i] = Byte.valueOf(temp[j]);
			site[index][i] = Integer.valueOf(temp[j+1]);
			jPos[index][i] = Integer.valueOf(temp[j+2]);
			jBinomP[index][i] = Double.valueOf(temp[j+3]);
			jSigSNPNumChr[index][i] = Integer.valueOf(temp[j+4]);
			likelyhood[index][i] = Double.valueOf(temp[j+5]);
			familyNum[index][i] = Byte.valueOf(temp[j+6]);
			familyCode[index][i] = Long.valueOf(temp[j+7]);
		}
	}

	private int getSize(String infileS) {
		int cnt = 0;
		try {
			BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
			String temp = br.readLine();
			temp = br.readLine();
			cnt++;
			String[] tem = temp.split("\\s+");
			long[] temTag = BaseEncoder.getLongArrayFromSeq(tem[0]);
			tagLengthInLong = temTag.length;
			while ((temp = br.readLine()) != null) cnt++;
		}
		catch (Exception e) {System.out.println(e.toString());}
		System.out.println(cnt + " tags int tagPGMap file " + infileS);
		return cnt;
	}

	public void writeTxtFile (String outfileS, boolean[] ifOut) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			bw.write(this.getHeader());
			bw.newLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				if (!ifOut[i]) continue;
				bw.write(this.getOutputString(i));
				bw.newLine();
			}
			bw.flush();
			bw.close();
			System.out.println("Txt tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	@Override
	public void writeTxtFile (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
			bw.write(this.getHeader());
			bw.newLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				bw.write(this.getOutputString(i));
				bw.newLine();
			}
			bw.flush();
			bw.close();
			System.out.println("Txt tagPGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	private String getOutputString (int index) {
		StringBuilder sb = new StringBuilder();	
		for (int i = 0; i < this.tagLengthInLong; i++) {
			sb.append(BaseEncoder.getSequenceFromLong(tags[i][index]));
		}
		sb.append("\t").append(String.valueOf(tagLength[index])).append("\t").append(String.valueOf(pChr[index])).append("\t").append(String.valueOf(pPos[index])).append("\t");
		sb.append(String.valueOf(ifMatch[index])).append("\t").append(String.valueOf(ifPerfectMatch[index])).append("\t").append(String.valueOf(ifUniqueMatch[index])).append("\t");
		sb.append(String.valueOf(this.tagTaxaCount[index])).append("\t").append(String.valueOf(this.jNum[index])).append("\t");
		for (int i = 0; i < jNum[index]; i++) {
			sb.append(String.valueOf(jChr[index][i])).append("\t");
			sb.append(String.valueOf(site[index][i])).append("\t");
			sb.append(String.valueOf(jPos[index][i])).append("\t");
			sb.append(String.valueOf(jBinomP[index][i])).append("\t");
			sb.append(String.valueOf(jSigSNPNumChr[index][i])).append("\t");
			sb.append(String.valueOf(likelyhood[index][i])).append("\t");
			sb.append(String.valueOf(familyNum[index][i])).append("\t");
			sb.append(String.valueOf(familyCode[index][i])).append("\t");
		}
		return sb.toString();
	}

	private String getHeader() {
		String header = "TestTag\tTagLength\tPChr\tPPos\tIfMatch\tIfPerfectMatch\tIfUniqueMatch\tTagTaxaCount\tFamilyGroupNum\t";
		for (int i = 0; i < familyGroupNum; i++) {
			String sufix = "-G"+String.valueOf(i+1);
			header = header + "LDChr" + sufix + "\t";
			header = header + "LDSite" + sufix + "\t";
			header = header + "LDPos" + sufix + "\t";
			header = header + "BinomP" + sufix + "\t";
			header = header + "ChrSig" + sufix + "\t";
			header = header + "Likelyhood" + sufix + "\t";
			header = header + "FamilyNum" + sufix + "\t";
			header = header + "FamilyCode" + sufix + "\t";
		}
		return header;
	}

    public void sortByPhysicalPosition () {
        GenericSorting.quickSort(0, this.getTagCount(), compByPhysicalPosition, this);
		System.out.println("TagPGMap is sorted by physical position");
    }
    
	public void sortByGeneticPosition () {
		GenericSorting.quickSort(0, this.getTagCount(), compByGeneticPosition, this);
		System.out.println("TagPGMap is sorted by genetic position");
	}

    IntComparator compByPhysicalPosition = new IntComparator() {
		public int compare(int a, int b) {
			if (pChr[a] == pChr[b]) {
				return pPos[a]-pPos[b];
			}
			else {
				return pChr[a]-pChr[b];
			}
		}
	};
    
	IntComparator compByGeneticPosition = new IntComparator() {
		public int compare(int a, int b) {
			if (jChr[a][0] == jChr[b][0]) {
				return jPos[a][0]-jPos[b][0];
			}
			else {
				return jChr[a][0]-jChr[b][0];
			}
		}
	};

	@Override
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
		int temInt;
		temInt = tagTaxaCount[index1]; tagTaxaCount[index1] = tagTaxaCount[index2]; tagTaxaCount[index2] = temInt;
        boolean temBoo;
        temBoo = ifMatch[index1]; ifMatch[index1] = ifMatch[index2]; ifMatch[index2] = temBoo;
        temBoo = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = temBoo;
        temBoo = ifUniqueMatch[index1]; ifUniqueMatch[index1] = ifUniqueMatch[index2]; ifUniqueMatch[index2] = temBoo;
		temByte = jNum[index1]; jNum[index1] = jNum[index2]; jNum[index2] = temByte;
		byte[] tempByte;
		tempByte = jChr[index1]; jChr[index1] = jChr[index2]; jChr[index2] = tempByte;
		int[] tempInt;
		tempInt = site[index1]; site[index1] = site[index2]; site[index2] = tempInt;
		tempInt = jPos[index1]; jPos[index1] = jPos[index2]; jPos[index2] = tempInt;
		double[] tempD;
		tempD = jBinomP[index1]; jBinomP[index1] = jBinomP[index2]; jBinomP[index2] = tempD;
		tempInt = jSigSNPNumChr[index1]; jSigSNPNumChr[index1] = jSigSNPNumChr[index2]; jSigSNPNumChr[index2] = tempInt;
		tempD = likelyhood[index1]; likelyhood[index1] = likelyhood[index2]; likelyhood[index2] = tempD;
		tempByte = familyNum[index1]; familyNum[index1] = familyNum[index2]; familyNum[index2] = tempByte;
		long[] tempLong;
		tempLong = familyCode[index1]; familyCode[index1] = familyCode[index2]; familyCode[index2] = tempLong;
		temByte = pChr[index1]; pChr[index1] = pChr[index2]; pChr[index2] = temByte;
		temInt = pPos[index1]; pPos[index1] = pPos[index2]; pPos[index2] = temInt;
		boolean tempB;
		tempB = ifMatch[index1]; ifMatch[index1] = ifMatch[index2]; ifMatch[index2] = tempB;
		tempB = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = tempB;
		tempB = ifUniqueMatch[index1]; ifUniqueMatch[index1] = ifUniqueMatch[index2]; ifUniqueMatch[index2] = tempB;
		
    }
}
