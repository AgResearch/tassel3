/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

/**
 *
 * @author Fei Lu
 */
public class VariableBinByTaxa {
    String[] taxaName;
    int[] tagTaxaCount;
    int[] chrID;
    int[][] binStart;
    int[][] binSize;
    short[][][] binCount;
    
    public VariableBinByTaxa () {
        
    }
    
    public VariableBinByTaxa (String genomeInfoFileS, String tagPMapFileS, String binSizeFileS, String tbtFileS) {
        this.creatBins(binSizeFileS, genomeInfoFileS);
        this.calBinCount(tagPMapFileS, tbtFileS);
        this.calTaxaCount();
    }
    
    public VariableBinByTaxa (String VariableBinByTaxaFileS) {
        this.readBinaryFile(VariableBinByTaxaFileS);
        this.calTaxaCount();
    }
    
    /*
     * select taxa from a table, made a new vbbt
     */
    public void selectTaxa (String taxaSelectionFileS) {
        Table t = new Table (taxaSelectionFileS);
        String[] candiTaxaName = new String[t.getRowNumber()];
        for (int i = 0; i < candiTaxaName.length; i++) {
            candiTaxaName[i] = t.content[i][0];
        }
        Arrays.sort(candiTaxaName);
        ArrayList<String> shareTaxaList = new ArrayList();
        for (int i = 0; i < this.getTaxaNum(); i++) {
            int hit = Arrays.binarySearch(candiTaxaName, this.taxaName[i]);
            if (hit < 0) continue;
            shareTaxaList.add(taxaName[i]);
        }
        String[] shareTaxa = shareTaxaList.toArray(new String[shareTaxaList.size()]);
        Arrays.sort(shareTaxa);
        int[] index = new int[shareTaxa.length];
        for (int i = 0; i < shareTaxa.length; i++) {
            for (int j = 0; j < taxaName.length; j++) {
                if (shareTaxa[i].equals(taxaName[j])) {
                    index[i] = j;
                    break;
                }
            }
        }
        short[][][] newBinCount = new short[chrID.length][][];
         for (int i = 0; i < newBinCount.length; i++) {
            newBinCount[i] = new short[binStart[i].length][shareTaxa.length];
            for (int j = 0; j < newBinCount[i].length; j++) {
                for (int k = 0; k < shareTaxa.length; k++) {
                    newBinCount[i][j][k] = binCount[i][j][index[k]];
                }
            }
        }
        taxaName = shareTaxa;
        binCount = newBinCount;
        this.calTaxaCount();
        
    }
    
    /*
     * merge samples of the same taxon
     */
    public void mergeByTaxaName () {
        this.lowerTaxaName();
        TreeSet<String> mergeNameSet = new TreeSet();
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            if (temp[0].equals("blank")) continue;
            if (temp[0].equals("Blank")) continue;
            if (temp[0].equals("unknown")) continue;
            mergeNameSet.add(temp[0]);
        }
        String[] uniqueName = mergeNameSet.toArray(new String[mergeNameSet.size()]);        
        Arrays.sort(uniqueName);
        ArrayList<Integer>[] referList = new ArrayList[uniqueName.length];
        String[] shortTaxaName = new String[taxaName.length];
        for (int i = 0; i < shortTaxaName.length; i++) {
            shortTaxaName[i] = taxaName[i].split(":")[0];
        }
        for (int i = 0; i < uniqueName.length; i++) {
            referList[i] = new ArrayList();  
            for (int j = 0; j < this.getTaxaNum(); j++) {
                if (uniqueName[i].equals(shortTaxaName[j])) {
                    referList[i].add(j);
                }
            }
        }
        short[][][] newBinCount = new short[chrID.length][][];
        for (int i = 0; i < newBinCount.length; i++) {
            newBinCount[i] = new short[binStart[i].length][uniqueName.length];
            for (int j = 0; j < newBinCount[i].length; j++) {
                for (int k = 0; k < uniqueName.length; k++) {
                    newBinCount[i][j][k] = 0;
                    for (int u = 0; u < referList[k].size(); u++) {
                        newBinCount[i][j][k] += binCount[i][j][referList[k].get(u)];
                    }
                }
            }
        }
        taxaName = uniqueName;
        binCount = newBinCount;
        this.calTaxaCount();
    }
    
    /*
     * add b73 ref, merge all of B73, B73(PI550473), B73Htrhm
     */
    public void addb73Taxon () {
        String[] B73 = new String[3];
        B73[0] = "B73"; B73[1] = "B73(PI550473)"; B73[2] = "B73Htrhm";
        ArrayList<Integer> indexList = new ArrayList();
        for (int i = 0; i < taxaName.length; i++) {
            String query = taxaName[i].split(":")[0];
            int hit = Arrays.binarySearch(B73, query);
            if (hit < 0) continue;
            indexList.add(i);
        }
        TreeSet<String> uniqueNameSet = new TreeSet();
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            if (temp[0].equals("blank")) continue;
            if (temp[0].equals("Blank")) continue;
            if (temp[0].equals("unknown")) continue;
            int hit = Arrays.binarySearch(B73, temp[0]);
            if (hit >= 0) continue;
            uniqueNameSet.add(taxaName[i]);
        }
        uniqueNameSet.add("b73");
        String[] uniqueName = uniqueNameSet.toArray(new String[uniqueNameSet.size()]);        
        Arrays.sort(uniqueName);
        ArrayList<Integer>[] referList = new ArrayList[uniqueName.length];
        for (int i = 0; i < uniqueName.length; i++) {
            if (uniqueName[i].equals("b73")) {
                referList[i] = indexList;
            }
            else {
                referList[i] = new ArrayList();  
                for (int j = 0; j < this.getTaxaNum(); j++) {
                    if (uniqueName[i].equals(taxaName[j])) {
                        referList[i].add(j);
                    }
                }
            }
        }
        short[][][] newBinCount = new short[chrID.length][][];
        for (int i = 0; i < newBinCount.length; i++) {
            newBinCount[i] = new short[binStart[i].length][uniqueName.length];
            for (int j = 0; j < newBinCount[i].length; j++) {
                for (int k = 0; k < uniqueName.length; k++) {
                    newBinCount[i][j][k] = 0;
                    for (int u = 0; u < referList[k].size(); u++) {
                        newBinCount[i][j][k] += binCount[i][j][referList[k].get(u)];
                    }
                }
            }
        }
        taxaName = uniqueName;
        binCount = newBinCount;
        this.calTaxaCount();
    }
    
    public double screePrintAverCount () {
        double averCount = 0;
        int binNum = 0;
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                binNum+=this.getTaxaNum();
            }
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < this.getTaxaNum(); k++) {
                    averCount += (double)binCount[i][j][k]/binNum;
                }
            }
        }
        System.out.println("Average count in each bin is " + averCount);
        return averCount;
    }
    
    private void lowerTaxaName () {
        for (int i = 0; i < this.getTaxaNum(); i++) {
            taxaName[i] = taxaName[i].toLowerCase();
        }
    }
    
    public int getTaxaNum () {
        return this.taxaName.length;
    }
    
    /*
     * order taxaName 
     */
    public void orderTaxa () {
        String[] newTaxaName = new String[taxaName.length];
        System.arraycopy(taxaName, 0, newTaxaName, 0, taxaName.length);
        Arrays.sort(newTaxaName);
        int[] index = new int[taxaName.length];
        for (int i = 0; i < taxaName.length; i++) {
            index[i] = Arrays.binarySearch(newTaxaName, taxaName[i]);
        }
        short[][][] newBinCount = new short[chrID.length][][];
        for (int i = 0; i < newBinCount.length; i++) {
            newBinCount[i] = new short[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    newBinCount[i][j][index[k]] = binCount[i][j][k];
                }
            }
        }
        taxaName = newTaxaName;
        binCount = newBinCount;
        this.calTaxaCount();
    }
    
    /*
     * find shared taxa name between key file and Janurary build
     */
    public void mkPaperTaxaNameTable (String paperKeyFileS, String paperTaxaNameFileS) {
        Table t = new Table (paperKeyFileS);
        String[] paperTaxaName = new String[t.getRowNumber()]; 
        for (int i = 0; i < t.getRowNumber(); i++) {
            String a = t.content[i][5];
            if (t.content[i][6].length() == 1) {
                a = a  + t.content[i][6];
            }
            else {
                a = a + t.content[i][6];
            }
            a = t.content[i][3]+":"+t.content[i][0]+":"+t.content[i][1]+":"+a;
            paperTaxaName[i] = a;
        }
        Arrays.sort(paperTaxaName);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (paperTaxaNameFileS), 65536);
            bw.write("TaxaName");
            bw.newLine();
            for (int i = 0; i < taxaName.length; i++) {
                int hit = Arrays.binarySearch(paperTaxaName, this.taxaName[i]);
                if (hit < 0) continue;
                if (taxaName[i].startsWith("blank")) continue;
                if (taxaName[i].startsWith("Blank")) continue;
                bw.write(this.taxaName[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /*
     * output taxa name to a file
     */
    public void mkTaxaNameTable (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("TaxaName");
            bw.newLine();
            for (int i = 0; i < taxaName.length; i++) {
                bw.write(taxaName[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void readBinaryFile (String infileS) {
        System.out.println("Reading "+infileS);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536*1024));
            int taxaNum = dis.readInt();
            taxaName = new String[taxaNum];
            for (int i = 0; i < taxaNum; i++) taxaName[i] = dis.readUTF();
            tagTaxaCount = new int[taxaNum];
            for (int i = 0; i < taxaNum; i++) tagTaxaCount[i] = dis.readInt();
            int chrNum = dis.readInt();
            chrID = new int[chrNum];
            for (int i = 0; i < chrNum; i++) chrID[i] = dis.readInt();
            int[] binNumOnChr = new int[chrNum];
            for (int i = 0; i < chrNum; i++) binNumOnChr[i] = dis.readInt();
            binStart = new int[chrNum][];
            binSize = new int[chrNum][];
            binCount = new short[chrNum][][];
            for (int i = 0; i < chrNum; i++) {
                binStart[i] = new int[binNumOnChr[i]];
                binSize[i] = new int[binNumOnChr[i]];
                binCount[i] = new short[binNumOnChr[i]][taxaNum];
                for (int j = 0; j < binNumOnChr[i]; j++) {
                    binStart[i][j] = dis.readInt();
                    binSize[i][j] = dis.readInt();
                }
            }
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < binNumOnChr[i]; j++) {
                    for (int k = 0; k < taxaNum; k++) {
                        binCount[i][j][k] = dis.readShort();
                    }
                }
            }
			dis.close();
            System.out.println("Binary BBT file read from " + infileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void writeBinaryFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536*1024));
            dos.writeInt(taxaName.length);
            for (int i = 0; i < taxaName.length; i++) dos.writeUTF(taxaName[i]);
            for (int i = 0; i < tagTaxaCount.length; i++) dos.writeInt(tagTaxaCount[i]);
            dos.writeInt(chrID.length);
            for (int i = 0; i < chrID.length; i++) {
                dos.writeInt(chrID[i]);
            }
            for (int i = 0; i < binStart.length; i++) dos.writeInt(binStart[i].length);
            for (int i = 0; i < binStart.length; i++) {
                for (int j = 0; j < binStart[i].length; j++) {
                    dos.writeInt(binStart[i][j]);
                    dos.writeInt(binSize[i][j]);
                }
            }
            for (int i = 0; i < binCount.length; i++) {
                for (int j = 0; j < binCount[i].length; j++) {
                    for (int k = 0; k < binCount[i][j].length; k++) {
                        dos.writeShort(binCount[i][j][k]);
                    }
                }
            }
            dos.flush();
            dos.close();
            System.out.println("Binary BBT file is written to " + outfileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /*
     * output txt format for debugging
     */
    public void writeTxtFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("Chr\tBinStart\tBinSize");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < binCount.length; i++) {
                for (int j = 0; j < binCount[i].length; j++) {
                    bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(binStart[i][j])+"\t"+String.valueOf(binSize[i][j]));
                    for (int k = 0; k < binCount[i][j].length; k++) {
                        bw.write("\t"+String.valueOf(binCount[i][j][k]));
                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            System.out.println("Txt BBT file is written to " + outfileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /*
     * calculate count in bins
     */
    private void calBinCount (String tagPMapFileS, String tbtFileS) {
        TagPMap pmap = new TagPMap (tagPMapFileS);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(tbtFileS), 65536*1024));
            int tagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
            this.iniMatrix(taxaNum);
			for (int i = 0; i < taxaNum; i++) {
				taxaName[i] = dis.readUTF();	
			}
            for (int i = 0; i < tagNum; i++) {
				long[] t = new long[2];
				for (int j = 0; j < tagLengthInLong; j++) {
					t[j] = dis.readLong();
				}
				byte temp = dis.readByte();
                if (pmap.getIfUniqueTag(t)) {
                    int[] chrPos = pmap.getBestChrPos(t);
                    int chrIndex = Arrays.binarySearch(chrID, chrPos[0]);
                    int binIndex = this.getBinIndexOnChr(chrIndex, chrPos[1]);
                    for (int j = 0; j < taxaNum; j++) {
                        binCount[chrIndex][binIndex][j] += dis.readByte() ;
                    }
                }
                else {
                    for (int j = 0; j < taxaNum; j++) {
                        temp = dis.readByte();
                    }
                }
				
				if (i%100000 ==0) System.out.println(String.valueOf(i) + " tags checked in TBT");
			}
			dis.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /*
     * calculate count in one taxon
     */
    private void calTaxaCount () {
        tagTaxaCount = new int[taxaName.length];
        for (int i = 0; i < tagTaxaCount.length; i++) tagTaxaCount[i] = 0;
        for (int i = 0; i < tagTaxaCount.length; i++) {
            for (int j = 0; j < binCount.length; j++) {
                for (int k = 0; k < binCount[j].length; k++) {
                    tagTaxaCount[i] += binCount[j][k][i];
                }
                    
            }
        }
    }
    
    private int getBinIndexOnChr (int chrIndex, int pos) {
        int hit = Arrays.binarySearch(binStart[chrIndex], pos);
        if (hit < 0) return -hit-2;
        return hit;
    }
    
    private void creatBins (String binSizeFileS, String genomeInfoFileS) {
        Table t = new Table (binSizeFileS);
        Table tg = new Table (genomeInfoFileS);
        chrID = new int[tg.getRowNumber()];
        int[] chrLength = new int[chrID.length];
        binStart = new int[chrID.length][];
        binSize = new int[chrID.length][];
        for (int i = 0; i < chrID.length; i++) {
            chrID[i] = Integer.valueOf(tg.content[i][0]);
            chrLength[i] = Integer.valueOf(tg.content[i][1]);
        }
        int[] count = new int[chrID.length];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int hit = Arrays.binarySearch(chrID, Integer.valueOf(t.content[i][0]));
            if (hit < 0) continue;
            count[hit]++;
        }
        for (int i = 0; i < count.length; i++) {
            binStart[i] = new int[count[i]];
            binSize[i] = new int[count[i]];
            count[i] = 0;
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int hit = Arrays.binarySearch(chrID, Integer.valueOf(t.content[i][0]));
            if (hit < 0) continue;
            binStart[hit][count[hit]] = Integer.valueOf(t.content[i][1]);
            count[hit]++;
        }
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart[i].length-1; j++) {
                binSize[i][j] = binStart[i][j+1]-binStart[i][j];
            }
            binSize[i][binSize[i].length-1] = chrLength[i] - binStart[i][binStart[i].length-1];
        }
        System.out.println("Bins are created");
    }
    
    /*
     * output bin size to based on count cutoff in bin
     */
    public void mkBinSizeFile (String tagCountFileS, String tagPMapFileS, String binSizeFileS, int taxaNum, double averageCount) {
        TagCounts tc = new TagCounts (tagCountFileS, FilePacking.Byte);
        TagPMap tmp = new TagPMap (tagPMapFileS);
        tmp.sortByPosition();
        int sum = 0;
        int currentChr = tmp.pChr[0][0];
        int uniqueTagCount = 0;
        ArrayList<String> binList = new ArrayList();
        if (tmp.pChr[0][0] == 1) {
            binList.add("1\t1");
        }
        for (int i = 0; i < tmp.getTagCount(); i++) {
            if (tmp.pChr[i][0] < 1) continue;
            if (tmp.pChr[i][0] > 10) continue;
            
            if (tmp.pChr[i][0] != currentChr) {
                sum = 0;
                binList.add(String.valueOf(tmp.pChr[i][0])+"\t1");
                currentChr = tmp.pChr[i][0];
            }
            long[] tag = tmp.getTag(i);
            if (!tmp.getIfUniqueTag(i)) continue;
            uniqueTagCount++;
            int index = tc.getTagIndex(tag);
            sum+=tc.getReadCount(index);
            double aver = (double)sum/taxaNum;
            if (aver < averageCount) continue;
            binList.add(String.valueOf(tmp.pChr[i][0])+"\t"+String.valueOf(tmp.pPos[i][0]+1));
            //System.out.println(aver);
            sum = 0;
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (binSizeFileS), 65536);
            bw.write("Chr\tBinStart");
            bw.newLine();
            for (int i = 0; i < binList.size(); i++) {
                bw.write(binList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(uniqueTagCount + " unique tags are used in read depth analysis");
    }
    
    private void iniMatrix (int taxaNum) {
        taxaName = new String[taxaNum];
        binCount = new short[binStart.length][][];
        for (int i = 0; i < binCount.length; i++) {
            binCount[i] = new short[binStart[i].length][taxaNum];
            for (int j = 0; j < binCount[i].length; j++) {
                for (int k = 0; k < taxaNum; k++) {
                    binCount[i][j][k] = 0;
                }
            }
        }
        tagTaxaCount = new int[taxaNum];
        for (int i = 0; i < taxaNum; i++) {
            tagTaxaCount[i] = 0;
        }
        System.out.println("Memory used: " + String.valueOf((Runtime.getRuntime().maxMemory()-Runtime.getRuntime().freeMemory())/1024/1024) + "M");
    }
}
