/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.colt.GenericSorting;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.genome.BaseEncoder;

/**
 * TagLenth info is not here
 * @author Fei Lu
 */
public class TagJMap extends AbstractTags {
	int familyGroupNum = 2;
	int[] tagTaxaCount;
	byte[] jNum;
	byte[][] jChr;
	int[][] site;
	int[][] jPos;
	double[][] jBinomP;
	int[][] jSigSNPNumChr; //significant SNPs on the chromosome
	double[][] likelyhood;
	byte[][] familyNum;
	long[][] familyCode;

	public TagJMap () {}

	public TagJMap (String mapFileS, boolean ifBinary) {
		if (ifBinary) this.readBinaryFile(mapFileS);
		else this.readTxtFile(mapFileS);
		this.sortByTag();
	}
	
	private void iniMatrix (int tagNum) {
		tags = new long[tagLengthInLong][tagNum];
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
		System.out.println("Matrix of tagGMap is initialized with " + tagNum + " tags");
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

	public void readBinaryFile (String mapFileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(mapFileS), 65536));
			int tagNum = dis.readInt();
			this.tagLengthInLong = dis.readInt();
			this.iniMatrix(tagNum);
			for (int i = 0; i < tagNum; i++) {
				for (int j = 0; j < tagLengthInLong; j++) {
					tags[j][i] = dis.readLong();
				}
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
			System.out.println("Binary tagGMapFile is read from " + mapFileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	public void readTxtFile (String mapFileS) {
		try {
			BufferedReader bw = new BufferedReader (new FileReader(mapFileS), 65536);
			String temp = bw.readLine();
			temp = bw.readLine();
			String[] tem = temp.split("\\s+");
			long[] temTag = BaseEncoder.getLongArrayFromSeq(tem[0]);
			tagLengthInLong = temTag.length;
			int n = 1;
			while ((temp = bw.readLine()) != null) n++;
			bw.close();
			this.iniMatrix(n);
			bw = new BufferedReader (new FileReader(mapFileS), 65536);
			bw.readLine();
			for (int i = 0; i < this.getTagCount(); i++) {
				temp = bw.readLine();
				this.readTxtRecord(temp, i);
			}
			bw.close();
			System.out.println("Txt tagGMapFile is read from " + mapFileS);
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
		tagTaxaCount[index] = Integer.valueOf(temp[1]);
		jNum[index] = Byte.valueOf(temp[2]);
		jChr[index] = new byte[jNum[index]];
		site[index] = new int[jNum[index]];
		jPos[index] = new int[jNum[index]];
		jBinomP[index] = new double[jNum[index]];
		jSigSNPNumChr[index] = new int[jNum[index]];
		likelyhood[index] = new double[jNum[index]];
		familyNum[index] = new byte[jNum[index]];
		familyCode[index] = new long[jNum[index]];
		for (int i = 0; i < jNum[index]; i++) {
			int j = i*8+3;
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

	public void writeBinaryFile (String outfileS) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
			dos.writeInt(this.getTagCount());
			dos.writeInt(this.tagLengthInLong);
			for (int i = 0; i < this.getTagCount(); i++) {
				for (int j = 0; j < this.tagLengthInLong; j++){
					dos.writeLong(tags[j][i]);
				}
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
			System.out.println("Binary tagGMap file is written in " + outfileS);
		}
		catch (Exception e) {
			System.err.println(e.toString());
      e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeTxtFile (String outfileS) {
		System.out.println("Not supported yet");
	}

	private String getHeader() {
		String header = "TestTag\tTagTaxaCount\tFamilyGroupNum\t";
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

	public void sortByTag() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
		System.out.println("TagGMap file is sorted by tags");
    }

	@Override
	public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp=tags[i][index1];
            tags[i][index1]=tags[i][index2];
            tags[i][index2]=temp;
        }
		int temInt;
		temInt = tagTaxaCount[index1]; tagTaxaCount[index1] = tagTaxaCount[index2]; tagTaxaCount[index2] = temInt;
		byte temByte;
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
    }
}
