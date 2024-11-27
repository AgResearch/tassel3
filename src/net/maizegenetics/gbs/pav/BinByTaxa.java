/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

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
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author Fei Lu
 */
public class BinByTaxa {   
    String[] taxaName;
    int[] tagTaxaCount;
    int[] chrID;
    int[][] binStart;
    short[][][] tagBinCount;
    
    
    public BinByTaxa (String genomeInfoFileS, String tagPMapFileS, String mergedTBTFileS) {
        this.creatBins(genomeInfoFileS, 50000);
        this.calBinCount(tagPMapFileS, mergedTBTFileS);
        this.calTaxaCount();
    }
    
    public BinByTaxa (String[] taxaName, int[] chrID, int[][] binStart, short[][][] tagBinCount) {
        this.taxaName = taxaName;
        this.chrID = chrID;
        this.binStart = binStart;
        this.tagBinCount = tagBinCount;
        this.calTaxaCount();
    }
    
    public BinByTaxa (String binByTaxaFileS) {
        if (binByTaxaFileS.endsWith("txt")) {
            this.readTxtFile(binByTaxaFileS);
        }
        else if (binByTaxaFileS.endsWith("bin")) {
            this.readBinaryFile(binByTaxaFileS);
        }
        else {
            System.out.println("Doesn't support this suffix/format, program will stop");
            System.exit(1);
        }
        this.calTaxaCount();
    }
    
    public int getTagCountAll () {
        int cnt = 0;
        for (int i = 0; i < tagTaxaCount.length; i++) {
            cnt+= this.tagTaxaCount[i];
        }
        return cnt;
    }
    
    public int getTagCountOfTaxa (int index) {
        return this.tagTaxaCount[index];
    }
    
    public int getChrNum () {
        return chrID.length;
    }
    
    public int getTaxaNum () {
        return taxaName.length;
    }
    
    private void renameB73 () {
        for (int i = 0; i < taxaName.length; i++) {
            if (taxaName[i].split(":")[0].equals("B73(PI550473)")) taxaName[i] = "B73:"+String.valueOf(i);
        }
    }
    
    public void selectTaxa (String taxaNameFileS) {
        Table t = new Table (taxaNameFileS);
        String[] newTaxaName = new String[t.getRowNumber()];
        for (int i = 0; i < newTaxaName.length; i++) {
            newTaxaName[i] = t.content[i][0];
        }
        Arrays.sort(newTaxaName);
        int[] index = new int[newTaxaName.length];
        for (int i = 0; i < newTaxaName.length; i++) {
            for (int j = 0; j < taxaName.length; j++) {
                if (newTaxaName[i].equals(taxaName[j])) {
                    index[i] = j;
                    break;
                }
            }
        }
        short[][][] newBinCount = new short[chrID.length][][];
         for (int i = 0; i < newBinCount.length; i++) {
            newBinCount[i] = new short[binStart[i].length][newTaxaName.length];
            for (int j = 0; j < newBinCount[i].length; j++) {
                for (int k = 0; k < newTaxaName.length; k++) {
                    newBinCount[i][j][k] = tagBinCount[i][j][index[k]];
                }
            }
        }
        taxaName = newTaxaName;
        tagBinCount = newBinCount;
        this.calTaxaCount();
    }
    
    public void mkTaxaCountFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("Taxa\tCount");
            bw.newLine();
            for (int i = 0; i < taxaName.length; i++) {
                bw.write(taxaName[i]+"\t"+String.valueOf(this.tagTaxaCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void lowerTaxaName () {
        for (int i = 0; i < taxaName.length; i++) taxaName[i] = taxaName[i].toLowerCase();
    }
    
    public void mergeBins (int binSize) {
        int currentSize = this.binStart[0][1] - this.binStart[0][0];
        if (!(binSize%currentSize == 0)) {
            System.out.println("binSize is not qualified for merging bins");
            System.exit(1);
        }
        int[][] newBinStart = new int[this.getChrNum()][];
        short[][][] newTagBinCount = new short[this.getChrNum()][][];
        for (int i = 0; i < binStart.length; i++) {
            int step = binSize/currentSize;
            int left = binStart[i].length % step;
            int n;
            if (left == 0) {
                n = binStart[i].length / step;
            }
            else {
                n = binStart[i].length / step + 1;
            }
            newBinStart[i] = new int[n];
            newTagBinCount[i] = new short[n][this.getTaxaNum()];
            if (left == 0) {
                for (int j = 0; j < newBinStart[i].length; j++) {
                    newBinStart[i][j] = binStart[i][step*j];
                }
            }
            else {
                for (int j = 0; j < newBinStart[i].length-1; j++) {
                    newBinStart[i][j] = binStart[i][step*j];
                }
                newBinStart[i][newBinStart[i].length-1] = binStart[i][binStart[i].length-left]; 
            }
        }
        for (int i = 0; i < tagBinCount.length; i++) {
            for (int j = 0; j < tagBinCount[i].length; j++) {
                int index = Arrays.binarySearch(newBinStart[i], binStart[i][j]);
                if (index < 0) index = -index-2;
                for (int k = 0; k < tagBinCount[i][j].length; k++) {
                    newTagBinCount[i][index][k] += tagBinCount[i][j][k];
                }
            }
        }
        binStart = newBinStart;
        tagBinCount = newTagBinCount;
    }
    
    public void addb73Taxa () {
        ArrayList<Integer> indexList = new ArrayList();
        for (int i = 0; i < taxaName.length; i++) {
            if (taxaName[i].split(":")[0].equals("B73(PI550473)")) indexList.add(i);
            if (taxaName[i].split(":")[0].equals("B73")) indexList.add(i);
        }
        
        TreeSet<String> uniqueNameSet = new TreeSet();
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            if (temp[0].equals("blank")) continue;
            if (temp[0].equals("Blank")) continue;
            if (temp[0].equals("unknown")) continue;
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
                        newBinCount[i][j][k] += tagBinCount[i][j][referList[k].get(u)];
                    }
                }
            }
        }
        taxaName = uniqueName;
        tagBinCount = newBinCount;
        this.calTaxaCount();
    }
    
    public void mergeTaxa () {
        TreeSet<String> uniqueNameSet = new TreeSet();
        this.renameB73();
        this.lowerTaxaName();
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            if (temp[0].equals("blank")) continue;
            if (temp[0].equals("Blank")) continue;
            if (temp[0].equals("unknown")) continue;
            uniqueNameSet.add(temp[0]);
        }
        String[] uniqueName = uniqueNameSet.toArray(new String[uniqueNameSet.size()]);        
        Arrays.sort(uniqueName);
        ArrayList<Integer>[] referList = new ArrayList[uniqueName.length];
        String[] nameHead = new String[taxaName.length];
        for (int i = 0; i < nameHead.length; i++) {
            nameHead[i] = taxaName[i].split(":")[0];
        }
        for (int i = 0; i < uniqueName.length; i++) {
            referList[i] = new ArrayList();
            for (int j = 0; j < nameHead.length; j++) {
                if (uniqueName[i].equals(nameHead[j])) {
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
                        newBinCount[i][j][k] += tagBinCount[i][j][referList[k].get(u)];
                    }
                }
            }
        }
        taxaName = uniqueName;
        tagBinCount = newBinCount;
        this.calTaxaCount();
    }
    
    public void readTxtFile (String infileS) {
        System.out.println("Reading " + infileS);
        try {
           BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
           String temp = br.readLine();
           String[] tem = temp.split("\t");
           int taxaNum = tem.length-2;
           taxaName = new String[taxaNum];
           for (int i = 0; i < taxaName.length; i++) taxaName[i] = tem[i+2];
           tagTaxaCount = new int[taxaNum];
           for (int i = 0; i < tagTaxaCount.length; i++) tagTaxaCount[i] = 0;
           int[] binNumOnChr = new int[100];
           for (int i = 0; i < binNumOnChr.length; i++) binNumOnChr[i] = 0;
           while ((temp = br.readLine()) != null) {
               temp = temp.substring(0, 10);
               tem = temp.split("\t");
               int index = Integer.valueOf(tem[0]) - 1;
               binNumOnChr[index]++;
           }
           int chrNum = 0;
           for (int i = 0; i < binNumOnChr.length; i++) {
               if (binNumOnChr[i] == 0) continue;
               chrNum++;
           }
           chrID = new int[chrNum];
           int[] binNum = new int[chrNum];
           int cnt = 0;
           for (int i = 0; i < binNumOnChr.length; i++) {
               if (binNumOnChr[i] == 0) continue;
               chrID[cnt] = i+1;
               binNum[cnt] = binNumOnChr[i];
               cnt++;
           }
           
           binStart = new int[chrNum][];
           tagBinCount = new short[chrNum][][];
           for (int i = 0; i < binStart.length; i++) {
               binStart[i] = new int[binNum[i]];
               tagBinCount[i] = new short[binNum[i]][taxaNum]; 
           }
           br = new BufferedReader (new FileReader(infileS), 65536);
           br.readLine();
           for (int i = 0; i < chrID.length; i++) {
               for (int j = 0; j < binStart[i].length; j++) {
                   temp = br.readLine();
                   tem = temp.split("\t");
                   binStart[i][j] = Integer.valueOf(tem[1]);
                   for (int k = 0; k < taxaNum; k++) {
                       tagBinCount[i][j][k] = Short.valueOf(tem[k+2]);
                   }
               }
           }
           br.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
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
            tagBinCount = new short[chrNum][][];
            for (int i = 0; i < chrNum; i++) {
                binStart[i] = new int[binNumOnChr[i]];
                tagBinCount[i] = new short[binNumOnChr[i]][taxaNum];
                for (int j = 0; j < binNumOnChr[i]; j++) binStart[i][j] = dis.readInt();
            }
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < binNumOnChr[i]; j++) {
                    for (int k = 0; k < taxaNum; k++) {
                        tagBinCount[i][j][k] = dis.readShort();
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
                }
            }
            for (int i = 0; i < tagBinCount.length; i++) {
                for (int j = 0; j < tagBinCount[i].length; j++) {
                    for (int k = 0; k < tagBinCount[i][j].length; k++) {
                        dos.writeShort(tagBinCount[i][j][k]);
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
    
    public void writeTxtFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("Chr\tBinStart");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < tagBinCount.length; i++) {
                for (int j = 0; j < tagBinCount[i].length; j++) {
                    bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(binStart[i][j]));
                    for (int k = 0; k < tagBinCount[i][j].length; k++) {
                        bw.write("\t"+String.valueOf(tagBinCount[i][j][k]));
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
    
    private void calTaxaCount () {
        tagTaxaCount = new int[taxaName.length];
        for (int i = 0; i < tagTaxaCount.length; i++) tagTaxaCount[i] = 0;
        for (int i = 0; i < tagTaxaCount.length; i++) {
            for (int j = 0; j < tagBinCount.length; j++) {
                for (int k = 0; k < tagBinCount[j].length; k++) {
                    tagTaxaCount[i] += tagBinCount[j][k][i];
                }
                    
            }
        }
    }
    
    private void calBinCount (String tagPMapFileS, String mergedTBTFileS) {
        TagPMap pmap = new TagPMap (tagPMapFileS);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(mergedTBTFileS), 65536*1024));
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
                        tagBinCount[chrIndex][binIndex][j] += dis.readByte() ;
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
    
    private int getBinIndexOnChr (int chrIndex, int pos) {
        int hit = Arrays.binarySearch(binStart[chrIndex], pos);
        if (hit < 0) return -hit-2;
        return hit;
    }
    
    private void iniMatrix (int taxaNum) {
        taxaName = new String[taxaNum];
        tagBinCount = new short[binStart.length][][];
        for (int i = 0; i < tagBinCount.length; i++) {
            tagBinCount[i] = new short[binStart[i].length][taxaNum];
            for (int j = 0; j < tagBinCount[i].length; j++) {
                for (int k = 0; k < taxaNum; k++) {
                    tagBinCount[i][j][k] = 0;
                }
            }
        }
        tagTaxaCount = new int[taxaNum];
        for (int i = 0; i < taxaNum; i++) {
            tagTaxaCount[i] = 0;
        }
        System.out.println("Memory used: " + String.valueOf((Runtime.getRuntime().maxMemory()-Runtime.getRuntime().freeMemory())/1024/1024) + "M");
    }
    
    private void creatBins (String genomeInfoFileS, int binSize) {
        Table t = new Table (genomeInfoFileS);
        chrID = new int[t.getRowNumber()];
        binStart = new int[t.getRowNumber()][];
        for (int i = 0; i < t.getRowNumber(); i++) {
            chrID[i] = Integer.valueOf(t.content[i][0]);
            int size = Integer.valueOf(t.content[i][1])/binSize + 1;
            binStart[i] = new int[size];
            for (int j = 0; j < size; j++) {
                binStart[i][j] = j*binSize + 1;
            }
        }
        System.out.println("Bins with size of " + binSize + "bp are created.");
        System.out.println("Memory used: " + String.valueOf((Runtime.getRuntime().maxMemory()-Runtime.getRuntime().freeMemory())/1024/1024) + "M");
    }

}
