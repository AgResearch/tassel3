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
public class RatioMatrix {
    String[] taxaName;
    String refName;
    int refIndex;
    double pavPTreshold;
    double cnvPTheshold;
    double occurence;
    int[] chrID;
    int[][] binStart;
    boolean[][] ifPAV;
    boolean[][] ifCNV;
    double[][][] ratioMatrix;
    double[][][] pMatrix;
        
    public RatioMatrix (BinByTaxa bbt, String refName, double fdr, double occurence) {
        this.iniMatrix(bbt, refName);
        this.calMatrix(bbt);
        this.transformMatrix();
        this.occurence = occurence;
        //this.discovery(pavP, cnvP, occurence);
        this.discovery(bbt, fdr, occurence);
    }
    
    public RatioMatrix (BinByTaxa bbt, String refName, double pavP, double cnvP, double occurence) {
        this.iniMatrix(bbt, refName);
        this.calMatrix(bbt);
        this.transformMatrix();
        this.occurence = occurence;
        this.discovery(pavP, cnvP, occurence);
    }
    
    public RatioMatrix (BinByTaxa bbt, String refName, String zscoreTableS, double zCutoff, double occurence) {
        this.iniMatrix(bbt, refName);
        this.calMatrix(bbt);
        this.transformMatrix();
        this.discoveryByZscore(zscoreTableS, zCutoff, occurence);
    }
    
    public RatioMatrix (BinByTaxa bbt, String refName, String pavPCutoffFileS, String cnvPCutoffFileS, double occurence, boolean ifPCut) {
        this.iniMatrix(bbt, refName);
        this.calMatrix(bbt);
        this.transformMatrix();
        this.occurence = occurence;
        if (ifPCut) {
            this.discoveryByPCut(pavPCutoffFileS, cnvPCutoffFileS, occurence);
        }
        else {
            this.discoveryByRCut(pavPCutoffFileS, cnvPCutoffFileS, occurence);
        }
    }
    
    public RatioMatrix (String CNVPMatrixFileS) {
        this.readBinaryFile(CNVPMatrixFileS);
    }
    
    public int getCNVLociAll () {
        int cnt = 0;
        for (int i = 0; i < this.ifCNV.length; i++) cnt+= this.getCNVLociOnChr(i);
        return cnt;
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
    
    public void discovery (BinByTaxa bbt, double fdr, double occurence) {
        int refBinCountCut = 200;
        int countTresh = (int)(taxaName.length * occurence);
        System.out.println("Discovering PAVs");
        double[] pThresh = this.getPAVPThreshOfTaxa(fdr);
        for (int i = 0; i < ifPAV.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < ifPAV[i].length; j++) {
                if (bbt.tagBinCount[i][j][refIndex] < refBinCountCut) continue;
                int cnt = 0;
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] < 0 && pMatrix[i][j][k] < pThresh[k]) cnt++;
                }
                if (cnt > countTresh) {
                    ifPAV[i][j] = true;
                }
                else {
                    ifPAV[i][j] = false;
                }
            }
       }
       pThresh = this.getCNVPThreshOfTaxa(fdr);
       for (int i = 0; i < ifCNV.length; i++) {
            for (int j = 0; j < ifCNV[i].length; j++) {
                int cnt = 0;
                if (bbt.tagBinCount[i][j][refIndex] < refBinCountCut) continue;
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > 0 && pMatrix[i][j][k] < pThresh[k]) cnt++;
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
    
    public void discoveryByRCut (String pavPCutoffFileS, String cnvPCutoffFileS, double occurence) {
        Table tpav = new Table (pavPCutoffFileS);
        Table tcnv = new Table (cnvPCutoffFileS);
        int countTresh = (int)(taxaName.length * occurence);
        System.out.println("Discovering PAVs and CNVs using p distribution along genome");
        int count = 0;
        for (int i = 0; i < ifPAV.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < ifPAV[i].length; j++) {
                int cnt = 0;
                double cut = Double.valueOf(tpav.content[count][2]);
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] < cut) cnt++;
                }
                if (cnt > countTresh) {
                    ifPAV[i][j] = true;
                }
                else {
                    ifPAV[i][j] = false;
                }
                count++;
            }
        }
        count = 0;
        for (int i = 0; i < ifCNV.length; i++) {
            for (int j = 0; j < ifCNV[i].length; j++) {
                int cnt = 0;
                double cut = Double.valueOf(tcnv.content[count][2]);
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > cut) cnt++;
                }
                if (cnt > countTresh) {
                    ifCNV[i][j] = true;
                }
                else {
                    ifCNV[i][j] = false;
                }
                count++;
            }
        }           
    }
    
    public void discoveryByZscore (String zscoreTableS, double zcutoff, double occurence) {
        Table t = new Table (zscoreTableS);
        double[] mean = new double[t.getRowNumber()];
        double[] se = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            mean[i] = Double.valueOf(t.content[i][3]);
            se[i] = Double.valueOf(t.content[i][2]);
        }
        int count = 0;
        int countTresh = (int)(taxaName.length * occurence);
        System.out.println("Discovering PAVs and CNVs using zscore along genome");
        for (int i = 0; i < ifPAV.length; i++) {
            System.out.println("Discoering PAVs and CNVs On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < ifPAV[i].length; j++) {
                int cntCNV = 0, cntPAV = 0;
                for (int k = 0; k < taxaName.length; k++) {
                    double z = (ratioMatrix[i][j][k]-mean[count])/se[i];
                    if (z > zcutoff) cntCNV++;
                    if (z < -zcutoff) cntPAV++;
                }
                if (cntCNV > countTresh) ifCNV[i][j] = true;
                if (cntPAV > countTresh) ifPAV[i][j] = true;
                count++;
            }
        }
    }
    
    public void discoveryByPCut (String pavPCutoffFileS, String cnvPCutoffFileS, double occurence) {
        Table tpav = new Table (pavPCutoffFileS);
        Table tcnv = new Table (cnvPCutoffFileS);
        int countTresh = (int)(taxaName.length * occurence);
        System.out.println("Discovering PAVs and CNVs using p distribution along genome");
        int count = 0;
        for (int i = 0; i < ifPAV.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < ifPAV[i].length; j++) {
                int cnt = 0;
                double cut = Double.valueOf(tpav.content[count][3]);
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] < 0 && pMatrix[i][j][k] < cut) cnt++;
                }
                if (cnt > countTresh) {
                    ifPAV[i][j] = true;
                }
                else {
                    ifPAV[i][j] = false;
                }
                count++;
            }
        }
        count = 0;
        for (int i = 0; i < ifCNV.length; i++) {
            for (int j = 0; j < ifCNV[i].length; j++) {
                int cnt = 0;
                double cut = Double.valueOf(tcnv.content[count][3]);
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > 0 && pMatrix[i][j][k] < cut) cnt++;
                }
                if (cnt > countTresh) {
                    ifCNV[i][j] = true;
                }
                else {
                    ifCNV[i][j] = false;
                }
                count++;
            }
        }           
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
    
    public void mkSVCountInRegion (int regionSize, String outfileS) {
        int currentSize = this.binStart[0][1] - this.binStart[0][0];
        if (!(regionSize%currentSize == 0)) {
            System.out.println("binSize is not qualified for merging bins");
            System.exit(1);
        }
        int[][] newBinStart = new int[this.getChrNum()][];
        int[][] pavCount = new int[this.getChrNum()][];
        int[][] cnvCount = new int[this.getChrNum()][];
        for (int i = 0; i < binStart.length; i++) {
            int step = regionSize/currentSize;
            int left = binStart[i].length % step;
            int n;
            if (left == 0) {
                n = binStart[i].length / step;
            }
            else {
                n = binStart[i].length / step + 1;
            }
            newBinStart[i] = new int[n];
            pavCount[i] = new int[n];
            cnvCount[i] = new int[n];
            
            if (left == 0) {
                for (int j = 0; j < newBinStart[i].length; j++) {
                    newBinStart[i][j] = binStart[i][step*j];
                    for (int k = step*j; k < step*j+step; k++) {
                        if (ifPAV[i][k]) pavCount[i][j]++;
                        if (ifCNV[i][k]) cnvCount[i][j]++;
                    }
                }
            }
            else {
                for (int j = 0; j < newBinStart[i].length-1; j++) {
                    newBinStart[i][j] = binStart[i][step*j];
                    for (int k = step*j; k < step*j+step; k++) {
                        if (ifPAV[i][k]) pavCount[i][j]++;
                        if (ifCNV[i][k]) cnvCount[i][j]++;
                    }
                }
                newBinStart[i][newBinStart[i].length-1] = binStart[i][binStart[i].length-left];
                for (int j = binStart[i].length-left; j < binStart[i].length; j++) {
                        if (ifPAV[i][j]) pavCount[i][newBinStart[i].length-1]++;
                        if (ifCNV[i][j]) cnvCount[i][newBinStart[i].length-1]++;
                }
            }
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Chr\tBinStart\tpavCount\tcnvCount");
            bw.newLine();
            for (int i = 0; i < newBinStart.length; i++) {
                for (int j = 0; j < newBinStart[i].length; j++) {
                    bw.write(String.valueOf(this.chrID[i])+"\t"+String.valueOf(newBinStart[i][j])+"\t"+String.valueOf(pavCount[i][j])+"\t"+String.valueOf(cnvCount[i][j]));
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
    
    public void mkCNVgenotypeByKmeans (String hapMapFileS) {
        boolean[][][] ifCNVgenotype = this.getIfCNVgenotypeByKmeans(new File(hapMapFileS).getParent());
        this.mkHapMapFile(hapMapFileS, ifCNVgenotype, false);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkCNVgenotypeByP (String hapMapFileS, double fdr) {
        boolean[][][] ifCNVgenotype = this.getIfCNVgenotypeByP(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifCNVgenotype, false);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkPAVgenotypeAverageRatioByP (String ratioFileS, double fdr) {
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByP(0, fdr);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(ratioFileS), 65536);
            bw.write("Chr\tBin\tDeletion");
            bw.newLine();
            for (int i = 0; i < this.ratioMatrix.length; i++) {
                for (int j = 0; j < this.ratioMatrix[i].length; j++) {
                    if (!ifPAV[i][j]) continue;
                    ArrayList<Double> list = new ArrayList();
                    for (int k = 0; k < this.getTaxaNum(); k++) {
                        if (!ifPAVgenotype[i][j][k]) continue;
                        list.add(ratioMatrix[i][j][k]);
                    }
                    Double[] ra = list.toArray(new Double[list.size()]);
                    double aver = 0;
                    for (int k = 0; k < ra.length; k++) aver = ra[k]/ra.length+aver;
                    bw.write(String.valueOf(this.chrID[i])+"\t"+String.valueOf(this.binStart[i][j])+"\t");
                    bw.write(String.valueOf(aver));
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
    
    public void mkPAVgenotypeMedianRatioByP (String ratioFileS, double fdr) {
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByP(0, fdr);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(ratioFileS), 65536);
            bw.write("Chr\tBin\tDeletion");
            bw.newLine();
            for (int i = 0; i < this.ratioMatrix.length; i++) {
                for (int j = 0; j < this.ratioMatrix[i].length; j++) {
                    if (!ifPAV[i][j]) continue;
                    ArrayList<Double> list = new ArrayList();
                    for (int k = 0; k < this.getTaxaNum(); k++) {
                        if (!ifPAVgenotype[i][j][k]) continue;
                        list.add(ratioMatrix[i][j][k]);
                    }
                    Double[] ra = list.toArray(new Double[list.size()]);
                    Arrays.sort(ra);
                    bw.write(String.valueOf(this.chrID[i])+"\t"+String.valueOf(this.binStart[i][j])+"\t");
                    bw.write(String.valueOf(ra[ra.length/2]));
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
    
    public void mkPAVgenotypeByP (String hapMapFileS, double fdr) {
        System.out.println(String.valueOf(this.getPAVLociAll())+" sites will be genotyped");
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByP(0, fdr);
        this.mkHapMapFile(hapMapFileS, ifPAVgenotype, true);
        System.out.println("PAV genotype is written to " + hapMapFileS);
    }
    
    public void mkPAVgenotypeByKmeans (String hapMapFileS) {
        System.out.println(String.valueOf(this.getPAVLociAll())+" sites will be genotyped");
        boolean[][][] ifPAVgenotype = this.getIfPAVgenotypeByKmeans(new File(hapMapFileS).getParent());
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
    
    private boolean[][][] getIfCNVgenotypeByKmeans (String dir) {
        boolean[][][] ifCNVgenotype = new boolean[this.chrID.length][][];
        String[] att = new String[2];
        att[0] = "r"; att[1] = "p";
        String wekaInstanceFileS = dir+"/in.arff";
        String outputFileS = dir+"/out.arff";
        for (int i = 0; i < this.getChrNum(); i++) {
            ifCNVgenotype[i] = new boolean[this.getBinNumOnChr(i)][this.getTaxaNum()];
            System.out.println("Genotyping at chromosome "+ String.valueOf(i+1));
            for (int j = 0; j < this.getBinNumOnChr(i); j++) {
                if (!ifCNV[i][j]) continue;
                double[][] value = new double[2][this.getTaxaNum()];
                for (int k = 0; k < this.getTaxaNum(); k++) {
                    value[0][k] = this.ratioMatrix[i][j][k];
                    value[1][k] = -Math.log(this.pMatrix[i][j][k]);
                    if (value[0][k] < 0) value[1][k] = -value[1][k];
                }
                WekaInstances wi = new WekaInstances(att, value);
                wi.writeNumericARFF(wekaInstanceFileS);
                try {
                    Runtime run = Runtime.getRuntime();
                    String cmd = "java weka.filters.unsupervised.attribute.AddCluster -W \"weka.clusterers.SimpleKMeans -A -N 2 -I 100\" -i " + wekaInstanceFileS + " -o " + outputFileS;                   
                    Process p = run.exec(cmd);
                    p.waitFor();
                }
                catch (Exception e) {
                    System.err.println(e.toString());
                    e.printStackTrace();
                    System.exit(1);
                }
                byte[] cluster = new WekaInstances().getCluster(outputFileS); 
                double max = 0;
                int index = 0;
                for (int k = 0; k < this.ratioMatrix[i][j].length; k++) {
                    if (this.ratioMatrix[i][j][k] > max) {
                        max = this.ratioMatrix[i][j][k];
                        index = k;
                    }
                }
                int cnvClusterValue;
                if (cluster[index] == 1) cnvClusterValue = 1;
                else cnvClusterValue = 2;
                
                for (int k = 0; k < this.ratioMatrix[i][j].length; k++) {
                    if (cluster[k] == cnvClusterValue) ifCNVgenotype[i][j][k] = true;
                    else ifCNVgenotype[i][j][k] = false;
                }
            }
        }
        System.out.println(this.getCNVLociAll());
        new File(wekaInstanceFileS).delete();
        new File(outputFileS).delete();
        return  ifCNVgenotype;
    }
    
    private boolean[][][] getIfPAVgenotypeByKmeans (String dir) {
        boolean[][][] ifPAVgenotype = new boolean[this.chrID.length][][];
        String[] att = new String[2];
        att[0] = "r"; att[1] = "p";
        String wekaInstanceFileS = dir+"/in.arff";
        String outputFileS = dir+"/out.arff";
        for (int i = 0; i < this.getChrNum(); i++) {
            ifPAVgenotype[i] = new boolean[this.getBinNumOnChr(i)][this.getTaxaNum()];
            System.out.println("Genotyping at chromosome "+ String.valueOf(i+1));
            for (int j = 0; j < this.getBinNumOnChr(i); j++) {
                if (!ifPAV[i][j]) continue;
                double[][] value = new double[2][this.getTaxaNum()];
                for (int k = 0; k < this.getTaxaNum(); k++) {
                    value[0][k] = this.ratioMatrix[i][j][k];
                    value[1][k] = 1-this.pMatrix[i][j][k];
                    if (value[0][k] < 0) value[1][k] = -value[1][k];
                }
                WekaInstances wi = new WekaInstances(att, value);
                wi.writeNumericARFF(wekaInstanceFileS);
                try {
                    Runtime run = Runtime.getRuntime();
                    String cmd = "java weka.filters.unsupervised.attribute.AddCluster -W \"weka.clusterers.SimpleKMeans -A -N 2 -I 100\" -i " + wekaInstanceFileS + " -o " + outputFileS;                   
                    Process p = run.exec(cmd);
                    p.waitFor();
                }
                catch (Exception e) {
                    System.err.println(e.toString());
                    e.printStackTrace();
                    System.exit(1);
                }
                byte[] cluster = new WekaInstances().getCluster(outputFileS); 
                double min = 0;
                int index = 0;
                for (int k = 0; k < this.ratioMatrix[i][j].length; k++) {
                    if (this.ratioMatrix[i][j][k] < min) {
                        min = this.ratioMatrix[i][j][k];
                        index = k;
                    }
                }
                int pavClusterValue;
                if (cluster[index] == 1) pavClusterValue = 1;
                else pavClusterValue = 2;
                
                for (int k = 0; k < this.ratioMatrix[i][j].length; k++) {
                    if (cluster[k] == pavClusterValue) ifPAVgenotype[i][j][k] = true;
                    else ifPAVgenotype[i][j][k] = false;
                }
            }
        }
        new File(wekaInstanceFileS).delete();
        new File(outputFileS).delete();
        return  ifPAVgenotype;
    }
    
    private boolean[][][] getIfCNVgenotypeByP (double ratioTresh, double fdr) {
        double[] pThresh = this.getCNVPThreshOfTaxa(fdr);
        boolean[][][] ifCNVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifCNVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifCNVgenotype[i].length; j++) {
                if (!ifCNV[i][j]) continue;
                for (int k = 0; k < ifCNVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] < ratioTresh) continue;
                    if (pMatrix[i][j][k] >= pThresh[k]) continue;
                    ifCNVgenotype[i][j][k] = true;
                }
            }
        }
        return ifCNVgenotype;
    }
    
    private boolean[][][] getIfPAVgenotypeByP (double ratioTresh, double fdr) {
        double[] pThresh = this.getPAVPThreshOfTaxa(fdr);
        boolean[][][] ifPAVgenotype = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifPAVgenotype[i] = new boolean[binStart[i].length][taxaName.length];
            for (int j = 0; j < ifPAVgenotype[i].length; j++) {
                if (!ifPAV[i][j]) continue;
                for (int k = 0; k < ifPAVgenotype[i][j].length; k++) {
                    if (ratioMatrix[i][j][k] > ratioTresh) continue;
                    if (pMatrix[i][j][k] >= pThresh[k]) continue;
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
    
    public void mkPvalueFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Chr\tBinStart");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < pMatrix.length; i++) {
                for (int j = 0; j < pMatrix[i].length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(binStart[i][j]));
                    for (int k = 0; k < pMatrix[i][j].length; k++) {
                        bw.write("\t"+String.valueOf(pMatrix[i][j][k])); 
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
    
    public void calMatrix (BinByTaxa bbt) {
        System.out.println("Calculating ratio matrix");
        for (int i = 0; i < pMatrix.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < pMatrix[i].length; j++) {
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    long[][] count = new long[2][2];
                    count[0][0] = bbt.tagBinCount[i][j][k];
                    count[0][1] = bbt.tagTaxaCount[k]-count[0][0];
                    count[1][0] = bbt.tagBinCount[i][j][refIndex];
                    count[1][1] = bbt.tagTaxaCount[refIndex]-count[1][0];
                    double nonRefRatio = (double)count[0][0]/bbt.tagTaxaCount[k];
                    double refRatio = (double)count[1][0]/bbt.tagTaxaCount[refIndex];
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
           
    private void iniMatrix (BinByTaxa bbt, String refName) {
        taxaName = bbt.taxaName;
        this.refName = refName;
        chrID = bbt.chrID;
        binStart = bbt.binStart;
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
        for (int i = 0; i < bbt.tagBinCount.length; i++) {
            pMatrix[i] = new double[bbt.tagBinCount[i].length][taxaName.length];
            ratioMatrix[i] = new double[bbt.tagBinCount[i].length][taxaName.length];
            for (int j = 0; j < bbt.tagBinCount[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    ratioMatrix[i][j][k] = -1;
                    pMatrix[i][j][k] = 0;
                }
            }
        }
    }
}
