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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math.stat.inference.TestUtils;

/**
 *
 * @author Fei Lu
 */
public class RatioMatrixVariableBin {
    String[] taxaName;
    String refName;
    int refIndex;
    double pavPTreshold;
    double cnvPTheshold;
    double occurence;
    int[] chrID;
    int[][] binStart;
    int[][] binSize;
    boolean[][] ifPAV;
    boolean[][] ifCNV;
    double[][][] ratioMatrix;
    double[][][] pMatrix;
    
    public RatioMatrixVariableBin (VariableBinByTaxa vbbt, String refName, double pavP, double cnvP, double occurence) {
        this.iniMatrix(vbbt, refName);
        this.calMatrix(vbbt);
        this.transformMatrix();
        this.occurence = occurence;
        this.discovery(pavP, cnvP, occurence);
    }
    
    public RatioMatrixVariableBin (String CNVPMatrixFileS) {
        this.readBinaryFile(CNVPMatrixFileS);
    }
    
    public int getCNVLociAll () {
        int cnt = 0;
        for (int i = 0; i < this.ifCNV.length; i++) cnt+= this.getCNVLociOnChr(i);
        return cnt;
    }
    
    public double getCNVPortion () {
        int len = 0;
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                if (ifCNV[i][j]) len += binSize[i][j];
            }
        }
        double portion = (double)len/this.getGenomeSize();
        System.out.println(portion + " is CNV");
        return portion;
    }
    
    public double getPAVPortion () {
        int len = 0;
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                if (ifPAV[i][j]) len += binSize[i][j];
            }
        }
        double portion = (double)len/this.getGenomeSize();
        System.out.println(portion + " is PAV");
        return portion;
    }
    
    public double getPAVOrCNVPortion () {
        int len = 0;
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                if (ifCNV[i][j] || ifPAV[i][j]) len += binSize[i][j];
            }
        }
        double portion = (double)len/this.getGenomeSize();
        System.out.println(portion + " is PAV or CNV");
        return portion;
    }
    
    public int getGenomeSize () {
        int len = 0;
        for (int i = 0; i < binSize.length; i++) {
            for (int j = 0; j < binSize[i].length; j++) {
                len += binSize[i][j];
            }
        }
        return len;
    }
    
    public int getCNVLociOnChr (int index) {
        int cnt = 0;
        for (int i = 0; i < binStart[index].length; i++) {
            if (ifCNV[index][i]) cnt++;
        }
        return cnt;
    }
    
    public int getPAVLociAll () {
        int cnt = 0;
        for (int i = 0; i < this.ifPAV.length; i++) cnt+= this.getPAVLociOnChr(i);
        return cnt;
    }
    
    public int getPAVAndCNVLociAll () {
        int cnt = 0;
        for (int i = 0; i < this.ifPAV.length; i++) {
            for (int j = 0; j < this.binStart[i].length; j++) {
                if (ifPAV[i][j] && ifCNV[i][j]) cnt++;
            }
        }
        return cnt;
    }
    
    public int getPAVOrCNVLociAll () {
        int cnt = 0;
        for (int i = 0; i < this.ifPAV.length; i++) {
            for (int j = 0; j < this.binStart[i].length; j++) {
                if (ifPAV[i][j] || ifCNV[i][j]) cnt++;
            }
        }
        return cnt;
    }
    
    public int getPAVLociOnChr (int index) {
        int cnt = 0;
        for (int i = 0; i < binStart[index].length; i++) {
            if (ifPAV[index][i]) cnt++;
        }
        return cnt;
    }
   
    public int getBinNumAll () {
        int cnt = 0;
        for (int i = 0; i < binStart.length; i++) cnt+= this.getBinNumOnChr(i);
        return cnt;
    }
    
    public int getBinNumOnChr (int index) {
        return binStart[index].length;
    }
    
    public int getChrNum () {
        return chrID.length;
    }
     
    public void discovery (double pavP, double cnvP, double occurence) {
        int countTresh = (int)(taxaName.length * occurence);
        System.out.println("Discovering PAVs");
        for (int i = 0; i < ifPAV.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < ifPAV[i].length; j++) {
                int cnt = 0;
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] < 0 && pMatrix[i][j][k] < pavP) cnt++;
                }
                if (cnt > countTresh) {
                    ifPAV[i][j] = true;
                }
                else {
                    ifPAV[i][j] = false;
                }
            }
       }
       for (int i = 0; i < ifCNV.length; i++) {
            for (int j = 0; j < ifCNV[i].length; j++) {
                int cnt = 0;
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > 0 && pMatrix[i][j][k] < cnvP) cnt++;
                }
                if (cnt > countTresh) {
                    ifCNV[i][j] = true;
                }
                else {
                    ifCNV[i][j] = false;
                }
            }
        }
    }
    
    public void mkCNVgenotypeByFDR (String hapMapFileS, double fdr) {
        boolean[][][] ifCNVgenotype = this.getIfCNVgenotypeByFDR(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifCNVgenotype, false);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkCNVgenotypeByP (String hapMapFileS, double fdr) {
        boolean[][][] ifCNVgenotype = this.getIfCNVgenotypeByP(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifCNVgenotype, false);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkPAVgenotypeByFDR (String hapMapFileS, double fdr) {
        System.out.println(String.valueOf(this.getPAVLociAll())+" sites will be genotyped");
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByFDR(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifPAVgenotype, true);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkPAVgenotypeByP (String hapMapFileS, double fdr) {
        System.out.println(String.valueOf(this.getPAVLociAll())+" sites will be genotyped");
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByP(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifPAVgenotype, true);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    private void mkHapMapFile (String hapMapFileS, boolean[][][] ifPositive, boolean pavOrcnv) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(hapMapFileS), 65536);
            StringBuilder sb = new StringBuilder ();
            sb.append("rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode");
            for (int i = 0; i < this.getTaxaNum(); i++) sb.append("\t").append(taxaName[i]);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < pMatrix.length; i++) {
                for (int j = 0; j < pMatrix[i].length; j++) {
                    if (pavOrcnv) {
                        if (!ifPAV[i][j]) continue;
                    }
                    else {
                        if (!ifCNV[i][j]) continue;
                    }
                    sb = new StringBuilder();
                    String ID = String.valueOf(i+1)+"-"+String.valueOf(j+1);
                    sb.append(ID).append("\tA/C\t").append(String.valueOf(i+1)).append("\t").append(String.valueOf(binStart[i][j]));
                    sb.append("\t+\tNA\tNA\tNA\tNA\tNA\tNA");
                    for (int k = 0; k < taxaName.length; k++) {
                        if (ifPositive[i][j][k]) {
                            sb.append("\t").append("C");
                        }
                        else {
                            sb.append("\t").append("A");
                        }
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
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
    
    private boolean[][][] getIfCNVgenotypeByFDR (double ratioTresh, double fdr) {
        double[] pThresh = this.getCNVPThreshOfTaxa(fdr);
        boolean[][][] ifCNVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifCNVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifCNVgenotype[i].length; j++) {
                if (!ifCNV[i][j]) continue;
                for (int k = 0; k < ifCNVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] < ratioTresh) continue;
                    if (pMatrix[i][j][k] >= pThresh[k]) continue;
                    //if (pMatrix[i][j][k] >= 0.001) continue;
                    ifCNVgenotype[i][j][k] = true;
                }
            }
        }
        return ifCNVgenotype;
    }
    
    private boolean[][][] getIfCNVgenotypeByP (double ratioTresh, double p) {
        
        boolean[][][] ifCNVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifCNVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifCNVgenotype[i].length; j++) {
                if (!ifCNV[i][j]) continue;
                for (int k = 0; k < ifCNVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] < ratioTresh) continue;
                    
                    if (pMatrix[i][j][k] >= 0.001) continue;
                    ifCNVgenotype[i][j][k] = true;
                }
            }
        }
        return ifCNVgenotype;
    }
    
    private boolean[][][] getIfPAVgenotypeByFDR (double ratioTresh, double fdr) {
        double[] pThresh = this.getPAVPThreshOfTaxa(fdr);
        boolean[][][] ifPAVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifPAVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifPAVgenotype[i].length; j++) {
                if (!ifPAV[i][j]) continue;
                for (int k = 0; k < ifPAVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] > ratioTresh) continue;
                    if (pMatrix[i][j][k] >= pThresh[k]) continue;
                    //if (pMatrix[i][j][k] >= 0.001) continue;
                    ifPAVgenotype[i][j][k] = true;
                }
            }
        }
        return ifPAVgenotype;
    }
    
    private boolean[][][] getIfPAVgenotypeByP (double ratioTresh, double p) {  
        boolean[][][] ifPAVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifPAVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifPAVgenotype[i].length; j++) {
                if (!ifPAV[i][j]) continue;
                for (int k = 0; k < ifPAVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] > ratioTresh) continue;              
                    if (pMatrix[i][j][k] >= 0.001) continue;
                    ifPAVgenotype[i][j][k] = true;
                }
            }
        }
        return ifPAVgenotype;
    }
    
    private double[] getPAVPThreshOfTaxa (double fdr) {
        double[] pThresh = new double[this.getTaxaNum()];
        for (int i = 0; i < this.getTaxaNum(); i++) {
            int cnt = 0;
            for (int j = 0; j < pMatrix.length; j++) {
                for (int k = 0; k < pMatrix[j].length; k++) {
                    if (ratioMatrix[j][k][i] < 0) cnt++;
                }
            }
            double[] pTaxa = new double[cnt];
            if (cnt ==0) {
                pThresh[i] = 1;
                continue;
            }
            cnt = 0;
            for (int j = 0; j < pMatrix.length; j++) {
                for (int k = 0; k < pMatrix[j].length; k++) {
                    if (ratioMatrix[j][k][i] < 0) {
                        pTaxa[cnt] = pMatrix[j][k][i];
                        cnt++;
                    }
                    
                }
            }
            Arrays.sort(pTaxa);
            pThresh[i] = pTaxa[(int)(pTaxa.length*fdr)];
        }
        return pThresh;
    }
    
    private double[] getCNVPThreshOfTaxa (double fdr) {
        double[] pThresh = new double[this.getTaxaNum()];
        for (int i = 0; i < this.getTaxaNum(); i++) {
            int cnt = 0;
            for (int j = 0; j < pMatrix.length; j++) {
                for (int k = 0; k < pMatrix[j].length; k++) {
                    if (ratioMatrix[j][k][i] > 0) cnt++;
                }
            }
            double[] pTaxa = new double[cnt];
            if (cnt ==0) {
                pThresh[i] = 1;
                continue;
            }
            cnt = 0;
            for (int j = 0; j < pMatrix.length; j++) {
                for (int k = 0; k < pMatrix[j].length; k++) {
                    if (ratioMatrix[j][k][i] > 0) {
                        pTaxa[cnt] = pMatrix[j][k][i];
                        cnt++;
                    }
                    
                }
            }
            Arrays.sort(pTaxa);
            pThresh[i] = pTaxa[(int)(pTaxa.length*fdr)];
        }
        return pThresh;
    }
    
    private double[] getPThreshOfTaxa (double fdr) {
        double[] pThresh = new double[this.getTaxaNum()];
        for (int i = 0; i < this.getTaxaNum(); i++) {
            double[] pTaxa = new double[this.getBinNumAll()];
            int cnt = 0;
            for (int j = 0; j < pMatrix.length; j++) {
                for (int k = 0; k < pMatrix[j].length; k++) {
                    pTaxa[cnt] = pMatrix[j][k][i];
                    cnt++;
                }
            }
            Arrays.sort(pTaxa);
            pThresh[i] = pTaxa[(int)(pTaxa.length*fdr)];
        }
        return pThresh;
    }
    
    public int getTaxaNum () {
        return taxaName.length;
    }
    
    public void writeBinaryFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536*1024));
            dos.writeInt(taxaName.length);
            dos.writeInt(chrID.length);
            for (int i = 0; i < binStart.length; i++) {
                dos.writeInt(binStart[i].length);
            }
            for (int i = 0; i < taxaName.length; i++) dos.writeUTF(taxaName[i]);
            dos.writeUTF(refName);
            dos.writeInt(refIndex);
            dos.writeDouble(this.pavPTreshold);
            dos.writeDouble(this.cnvPTheshold);
            dos.writeDouble(this.occurence);
            for (int i = 0; i < chrID.length; i++) dos.writeInt(chrID[i]);
            for (int i = 0; i < binStart.length; i++) {
                for (int j = 0; j < binStart[i].length; j++) {
                    dos.writeInt(binStart[i][j]);
                    dos.writeBoolean(ifPAV[i][j]);
                    dos.writeBoolean(ifCNV[i][j]);
                }
            }
            for (int i = 0; i < pMatrix.length; i++) {
                for (int j = 0; j < pMatrix[i].length; j++) {
                    for (int k = 0; k < pMatrix[i][j].length; k++) {
                        dos.writeDouble(ratioMatrix[i][j][k]);
                        dos.writeDouble(pMatrix[i][j][k]);
                    }
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void readBinaryFile (String infileS) {
        try {
             DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536*1024));
             int taxaNum = dis.readInt();
             int chrNum = dis.readInt();
             int[] binNum = new int[chrNum];
             for (int i = 0; i < binNum.length; i++) binNum[i] = dis.readInt();
             this.iniMatrix(taxaNum, chrNum, binNum);
             for (int i = 0; i < taxaName.length; i++) taxaName[i] = dis.readUTF();
             refName = dis.readUTF();
             refIndex = dis.readInt();
             this.pavPTreshold = dis.readDouble();
             this.cnvPTheshold = dis.readDouble();
             this.occurence = dis.readDouble();
             for (int i = 0; i < chrID.length; i++) chrID[i] = dis.readInt();
             for (int i = 0; i < binStart.length; i++) {
                 for (int j = 0; j < binStart[i].length; j++) {
                     binStart[i][j] = dis.readInt();
                     ifPAV[i][j] = dis.readBoolean();
                     ifCNV[i][j] = dis.readBoolean();
                 }
             }
             for (int i = 0; i < pMatrix.length; i++) {
                 for (int j = 0; j < pMatrix[i].length; j++) {
                     for (int k = 0; k < pMatrix[i][j].length; k++) {
                         ratioMatrix[i][j][k] = dis.readDouble();
                         pMatrix[i][j][k] = dis.readDouble();
                     }
                 }
             }
        }
        catch(Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void iniMatrix (int taxaNum, int chrNum, int[] binNum) {
        taxaName = new String[taxaNum];
        chrID = new int[chrNum];
        binStart = new int[chrNum][];
        ifPAV = new boolean[chrNum][];
        ifCNV = new boolean[chrNum][];
        ratioMatrix = new double[chrNum][][];
        pMatrix = new double[chrNum][][];
        for (int i = 0; i < binStart.length; i++) {
            binStart[i] = new int[binNum[i]];
            ifPAV[i] = new boolean[binNum[i]];
            ifCNV[i] = new boolean[binNum[i]];
            ratioMatrix[i] = new double[binNum[i]][taxaNum];
            pMatrix[i] = new double[binNum[i]][taxaNum];
        }
    }
    
    public void transformMatrix () {
        System.out.println("Transforming ratio matrix");
        for (int i = 0; i <  pMatrix.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < pMatrix[i].length; j++) {
                double min = 0;;
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] < 0) {
                        ratioMatrix[i][j][k] = 0;
                    }
                    else if (ratioMatrix[i][j][k] > 0)  {
                        ratioMatrix[i][j][k] = Math.log(ratioMatrix[i][j][k])/Math.log(2);
                        if (ratioMatrix[i][j][k] < min) min = ratioMatrix[i][j][k];
                    }
                    if (Double.isNaN(pMatrix[i][j][k])) pMatrix[i][j][k] = 1;
                    if (pMatrix[i][j][k] == 0) pMatrix[i][j][k] = Double.MIN_VALUE;
                }
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] == 0) {
                        if (pMatrix[i][j][k] < 0.05) {
                            ratioMatrix[i][j][k] = min;
                        }
                        else {
                            ratioMatrix[i][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    
    public void displayFDRpCutoff (double fdr) {
        this.getPAVCutoff(fdr);
        this.getCNVCutoff(fdr);
    }
    
    public double getPAVCutoff (double fdr) {
        double[] pArray = new double[5000];
        for (int i = 0; i < pArray.length; i++) {
             int j = (int)(Math.random()*this.getChrNum());
             int k = (int)(Math.random()*this.pMatrix[j].length);
             int u = (int)(Math.random()*this.getTaxaNum());
             if (ratioMatrix[j][k][u] >= 0) {
                 i--;
                 continue;
             }
             pArray[i] = pMatrix[j][k][u];
        }
        Arrays.sort(pArray);
        int index = (int)Math.floor((double)(fdr) * pArray.length);
        System.out.println("When FDR = " + String.valueOf(fdr) + ", PAV cutoff is " + String.valueOf(pArray[index]));
        return pArray[index];
    }
    
    public double getCNVCutoff (double fdr) {
        double[] pArray = new double[5000];
        for (int i = 0; i < pArray.length; i++) {
             int j = (int)(Math.random()*this.getChrNum());
             int k = (int)(Math.random()*this.pMatrix[j].length);
             int u = (int)(Math.random()*this.getTaxaNum());
             if (ratioMatrix[j][k][u] >= 0) {
                 i--;
                 continue;
             }
             pArray[i] = pMatrix[j][k][u];
        }
        Arrays.sort(pArray);
        int index = (int)Math.floor((double)(fdr) * pArray.length);
        System.out.println("When FDR = " + String.valueOf(fdr) + ", CNV cutoff is " + String.valueOf(pArray[index]));
        return pArray[index];
    }
    
    public double[] getPvalue (boolean ifPAV) {
        int limit = 20000;
        double[] pArray = new double[limit];
        for (int i = 0; i < limit; i++) {
            int j = (int)Math.floor(pMatrix.length*Math.random());
            int k = (int)Math.floor(pMatrix[j].length*Math.random());
            int u = (int)Math.floor(pMatrix[j][k].length*Math.random());
            if (pMatrix[j][k][u] == 0) {
                i--;
                continue;
            }
            if (ifPAV) {
                if (pMatrix[j][k][u] > 0) {
                    i--;
                    continue;
                }
            }
            else {
                if (pMatrix[j][k][u] < 0) {
                    i--;
                    continue;
                }
            }
            pArray[i] = pMatrix[j][k][u];
        }
        return pArray;
    }
    
    public void mkMedianRatioFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Chr\tBin\tDeletion\tCNV");
            bw.newLine();
            for (int i = 0; i < ratioMatrix.length; i++) {
                for (int j = 0; j < ratioMatrix[i].length; j++) {
                    ArrayList<Double> dList = new ArrayList();
                    ArrayList<Double> cList = new ArrayList();
                    for (int k = 0; k < ratioMatrix[i][j].length; k++) {
                        if (ratioMatrix[i][j][k] > 0) {
                            cList.add(ratioMatrix[i][j][k]);
                        }
                        else if (ratioMatrix[i][j][k] < 0) {
                            dList.add(ratioMatrix[i][j][k]);
                        }
                    }
                    Double[] d = dList.toArray(new Double[dList.size()]);
                    Double[] c = cList.toArray(new Double[cList.size()]);
                    Arrays.sort(d);
                    Arrays.sort(c);
                    bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(this.binStart[i][j])+"\t");
                    if (d.length == 0) bw.write("0\t");
                    else bw.write(String.valueOf(d[d.length/2])+"\t");
                    if (c.length == 0) bw.write("0");
                    else bw.write(String.valueOf(c[c.length/2]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkRatioFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Chr\tBinStart");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < ratioMatrix.length; i++) {
                for (int j = 0; j < ratioMatrix[i].length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(binStart[i][j]));
                    for (int k = 0; k < ratioMatrix[i][j].length; k++) {
                        bw.write("\t"+String.valueOf(ratioMatrix[i][j][k])); 
                    }
                    bw.newLine();
                }
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
    
    public void writeTxtFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Chr\tBinStart");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]+"\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < pMatrix.length; i++) {
                for (int j = 0; j < pMatrix[i].length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(binStart[i][j]));
                    for (int k = 0; k < pMatrix[i][j].length; k++) {
                        bw.write("\t"+String.valueOf(ratioMatrix[i][j][k])+"\t"+String.valueOf(pMatrix[i][j][k])); 
                    }
                    bw.newLine();
                }
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
    
    public void calMatrix (VariableBinByTaxa vbbt) {
        System.out.println("Calculating ratio matrix");
        for (int i = 0; i < pMatrix.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < pMatrix[i].length; j++) {
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    long[][] count = new long[2][2];
                    count[0][0] = vbbt.binCount[i][j][k];
                    count[0][1] = vbbt.tagTaxaCount[k]-count[0][0];
                    count[1][0] = vbbt.binCount[i][j][refIndex];
                    count[1][1] = vbbt.tagTaxaCount[refIndex]-count[1][0];
                    double nonRefRatio = (double)count[0][0]/vbbt.tagTaxaCount[k];
                    double refRatio = (double)count[1][0]/vbbt.tagTaxaCount[refIndex];
                    if (refRatio > 0) ratioMatrix[i][j][k] = nonRefRatio/refRatio;
                    else {
                        ratioMatrix[i][j][k] = -1;
                    }
                    try {
                        pMatrix[i][j][k] = TestUtils.chiSquareTest(count);
                    }
                    catch (Exception e) {
                        System.err.println(e.toString());
                        e.printStackTrace();
                        System.exit(1);
                    }
                }
            }
        }
    }
           
    private void iniMatrix (VariableBinByTaxa vbbt, String refName) {
        taxaName = vbbt.taxaName;
        this.refName = refName;
        chrID = vbbt.chrID;
        binStart = vbbt.binStart;
        binSize = vbbt.binSize;
        refIndex = Arrays.binarySearch(taxaName, refName);
        ratioMatrix = new double[chrID.length][][];
        pMatrix = new double[chrID.length][][];
        ifPAV = new boolean[binStart.length][];
        ifCNV = new boolean[binStart.length][];
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart.length; j++) {
                ifPAV[i] = new boolean[binStart[i].length];
                ifCNV[i] = new boolean[binStart[i].length];
            }
        }
        for (int i = 0; i < vbbt.binCount.length; i++) {
            pMatrix[i] = new double[vbbt.binCount[i].length][taxaName.length];
            ratioMatrix[i] = new double[vbbt.binCount[i].length][taxaName.length];
            for (int j = 0; j < vbbt.binCount[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    ratioMatrix[i][j][k] = -1;
                    pMatrix[i][j][k] = 0;
                }
            }
        }
    }
}
