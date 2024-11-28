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
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.genome.BaseEncoder;

/**
 * Chromosome from 0 to 12, 0 is unanchored, 11 is Mt, 12 is Pt, sChr and eChr are essential chromosome range of maize
 * Tags don't align to ref are excluded
 * @author Fei Lu
 */
public class TagPMap extends AbstractTags {
    byte[] hitNum;
    byte[][] pChr;
    int[][] pPos;
    boolean[][] ifPerfectMatch;
    static byte sChr = 1;
    static byte eChr = 10;
    
    public TagPMap (AbstractTags tc, String samFileS) {
        this.tagLengthInLong = tc.getTagSizeInLong();
        this.iniMatrix(this.getSize(samFileS));
        this.readSam(samFileS, tc);
        this.sortByTag();
    }
    
    public TagPMap (String tagPMapFileS) {
        this.readTagPMap(tagPMapFileS);
    }
    
    public int getPerfectHitNum (int index) {
        int cnt = 0;
        for (int i = 0; i < ifPerfectMatch[index].length; i++) {
            if (ifPerfectMatch[index][i]) cnt++;
        }
        return cnt;
    }
    
    public int getHitNum (int index) {
        return this.hitNum[index];
    }
    
    public int getHitNum (long[] tag) {
        int index = this.getTagIndex(tag);
        return this.getHitNum(index);
    }
    
    public byte getBestChr (int index) {
        return this.pChr[index][0];
    }
    
    public int getBestPos (int index) {
        return this.pPos[index][0];
    }
    
    public int[] getBestChrPos (int index) {
        if (hitNum[index] == 0) return null;
        int[] a = new int[2];
        a[0] = pChr[index][0];
        a[1] = pPos[index][0];
        return a;
    }
    
    public int[] getBestChrPos (long[] tag) {
        int index = this.getTagIndex(tag);
        if (index < 0) return null;
        return this.getBestChrPos(index);
    }
    
    public boolean getIfRefUniqueSingleBestTag (int index) {
        if (this.hitNum[index] == 1) {
            if (ifPerfectMatch[index][0] == true) return true;
        }
        else {
            if (ifPerfectMatch[index][0] == true && ifPerfectMatch[index][1] == false) return true;
        }
        return false;
    }
    
    public boolean getIfRefUniqueSingleBestTag (long[] tag) {
        int index = this.getTagIndex(tag);
        if (index < 0) return false;
        return getIfRefUniqueSingleBestTag(index);
    }
    
    public boolean getIfRefUniqueTag (int index) {
        if (this.hitNum[index] == 1 && ifPerfectMatch[index][0] == true) return true;
        return false;
    }
    
    public boolean getIfRefUniqueTag (long[] tag) { //false also means the tested tag is not included in pMap
        int index = this.getTagIndex(tag);
        if (index < 0) return false;
        return getIfRefUniqueTag(index);
    }
    
    public boolean getIfUniqueTag (int index) {
        if (this.hitNum[index] == 1 && pChr[index][0] >= TagPMap.sChr && pChr[index][0] <= TagPMap.eChr) return true;
        return false;
    }
    
    public boolean getIfUniqueTag (long[] tag) {//false also means the tested tag is not included in pMap
        int index = this.getTagIndex(tag);
        if (index < 0) return false;
        return this.getIfUniqueTag(index);
    }
    
    public void readTagPMap (String tagPMapFileS) {
        try {
            DataInputStream dis = new DataInputStream (new BufferedInputStream (new FileInputStream(tagPMapFileS), 65536));
            int tagNum = dis.readInt();
            this.tagLengthInLong = dis.readInt();
            this.iniMatrix(tagNum);
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    tags[j][i] = dis.readLong();
                }
                tagLength[i] = dis.readByte();
                hitNum[i] = dis.readByte();
                this.iniSubMatrix(hitNum[i], i);
                for (int j = 0; j < hitNum[i]; j++) {
                    pChr[i][j] = dis.readByte();
                    pPos[i][j] = dis.readInt();
                    ifPerfectMatch[i][j] = dis.readBoolean();
                }
            }
            dis.close();
            System.out.println("Binary tagPMapFile is read from " + tagPMapFileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void writeTagPMap (String tagPMapFileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tagPMapFileS), 65536));
            dos.writeInt(this.getTagCount());
            dos.writeInt(this.tagLengthInLong);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++){
                    dos.writeLong(tags[j][i]);
                }
                dos.writeByte(this.tagLength[i]);
                dos.writeByte(this.hitNum[i]);
                for (int j = 0; j < this.hitNum[i]; j++) {
                    dos.writeByte(pChr[i][j]);
                    dos.writeInt(pPos[i][j]);
                    dos.writeBoolean(ifPerfectMatch[i][j]);
                }
            }
            dos.flush();
            dos.close();
            System.out.println("Binary tagPMapFile is written in " + tagPMapFileS);
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private int getSize(String samFileS) {
        TreeSet<Integer> idSet = new TreeSet();
        try {
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                String[] tem = temp.split("\t");
                if (tem[2].startsWith("*")) continue;
                idSet.add(Integer.valueOf(tem[0]));
            }
        }
        catch (Exception e) {System.out.println(e.toString());}
        System.out.println(idSet.size() + " tags aligned to reference genome in sam file");
        return idSet.size();
    }
    
    private void readSam (String samFileS, AbstractTags tc) {
        System.out.println("Start reading sam file");
        String temp = null;
        try {
            int tagIndex = 0;
            BufferedReader br = new BufferedReader (new FileReader(samFileS), 65536);
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
                        tagIndex = this.importSAMRecord(tagAlnList, tagIndex, tagIDList.get(0), tc);
                        
                        tagAlnList.clear();
                        tagIDList.clear();
                        tagAlnList.add(temp);
                        tagIDList.add(Integer.valueOf(tem[0]));
                    }
                }
            }
            tagIndex = this.importSAMRecord(tagAlnList, tagIndex, tagIDList.get(0), tc);
            System.out.println("Imported physical position of " + String.valueOf(tagIndex) + " tags");
        }
        catch (Exception e) {System.out.println(e.toString()+"\t"+temp);}
        System.out.println("Sam file read. Physical position of tags imported");
    }
    
    private int importSAMRecord (ArrayList<String> tagAlnList, int tagIndex, int tagID, AbstractTags tc) {
        String[] temp = tagAlnList.get(0).split("\t");
        if (temp[2].startsWith("*")) return tagIndex;
        byte chr = Byte.valueOf(temp[2]);
        //if (chr < this.sChr || chr > this.eChr) return tagIndex;
        
        long[] tempTag = tc.getTag(tagID-1);
        for (int i = 0; i < tagLengthInLong; i++) {
            tags[i][tagIndex] = tempTag[i];
        }
        tagLength[tagIndex] = (byte)tc.getTagLength(tagID-1);
        
        String[] tagAln = tagAlnList.toArray(new String[tagAlnList.size()]);
        hitNum[tagIndex] = (byte)tagAln.length;
        pChr[tagIndex] = new byte[tagAln.length];
        pPos[tagIndex] = new int[tagAln.length];
        ifPerfectMatch[tagIndex] = new boolean[tagAln.length];
        for (int i = 0; i < tagAln.length; i++) {
            temp = tagAln[i].split("\t");
            pChr[tagIndex][i] = Byte.valueOf(temp[2]);
            pPos[tagIndex][i] = Integer.valueOf(temp[3]);
            if (temp[5].equals(String.valueOf(tc.getTagLength(tagID-1))+"M") && tagAln[i].contains("NM:i:0")) {
                ifPerfectMatch[tagIndex][i] = true;
            }
            else {
                ifPerfectMatch[tagIndex][i] = false;
            }
        }
        return ++tagIndex;
    }
    
    @Override
    public void swap(int index1, int index2) {
        super.swap(index1, index2);
        byte tempHitNum;
        tempHitNum = hitNum[index1]; hitNum[index1] = hitNum[index2]; hitNum[index2] = tempHitNum;
        byte[] tempChr;
        tempChr = pChr[index1]; pChr[index1] = pChr[index2]; pChr[index2] = tempChr;
        int[] tempPos;
        tempPos = pPos[index1]; pPos[index1] = pPos[index2]; pPos[index2] = tempPos;
        boolean[] tempIfPerfectMatch;
        tempIfPerfectMatch = ifPerfectMatch[index1]; ifPerfectMatch[index1] = ifPerfectMatch[index2]; ifPerfectMatch[index2] = tempIfPerfectMatch;
    }
    
    public void sortByTag() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
        System.out.println("TagPMap file is sorted by tags");
    }
    
    public void sortByPosition () {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), compByPosition, this);
        System.out.println("Position index sort end.");
        System.out.println("TagPMap file is sorted by position");
    }
    
    private void iniMatrix (int tagNum) {
        tagLengthInLong = 2;
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        hitNum = new byte[tagNum];
        pChr = new byte[tagNum][];
        pPos = new int[tagNum][];
        ifPerfectMatch = new boolean[tagNum][];
        System.out.println("Matrix of tagPMap is initialized with " + tagNum + " tags");
    }
    
    private void iniSubMatrix (int hitCount, int index) {
        pChr[index] = new byte[hitCount];
        pPos[index] = new int[hitCount];
        ifPerfectMatch[index] = new boolean[hitCount];
    }
    
    IntComparator compByPosition = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            if (pChr[a][0] < pChr[b][0]) {
                return -1;
            }
            else if (pChr[a][0] > pChr[b][0]) {
                return 1;
            }
            else {
                return pPos[a][0]-pPos[b][0];
            }
        }
    };
}
