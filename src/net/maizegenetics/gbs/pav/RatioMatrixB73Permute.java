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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TestUtils;

/**
 *
 * @author Fei Lu
 */
public class RatioMatrixB73Permute {
    String[] taxaName;
    String refName;
    int refIndex;
    int[] chrID;
    int[][] binStart;
    double[][][] ratioMatrix;
    double[][][] pMatrix;
           
    public RatioMatrixB73Permute (BinByTaxa bbt, String refName) {
        this.iniMatrix(bbt, refName);
        this.calMatrix(bbt);
        this.transformMatrix();
    }
    
    public void mkGenotype (String pavCutoffFileS, String cnvCutoffFileS, String pavGenotypeFileS, String cnvGenotypeFileS, double maf) {
        boolean[][][] ifPAV = this.discoverPAV(pavCutoffFileS);
        boolean[][][] ifCNV = this.discoverCNV(cnvCutoffFileS);
        this.mkHapMapFile(pavGenotypeFileS, ifPAV, maf);
        this.mkHapMapFile(cnvGenotypeFileS, ifCNV, maf);
    }
    
    public void mkGenotype (double pavCutoff, double cnvCutoff, String pavGenotypeFileS, String cnvGenotypeFileS, double maf) {
        boolean[][][] ifPAV = this.discoverPAV(pavCutoff);
        boolean[][][] ifCNV = this.discoverCNV(cnvCutoff);
        String[] pavId = this.mkHapMapFile(pavGenotypeFileS, ifPAV, maf);
        String[] cnvId = this.mkHapMapFile(cnvGenotypeFileS, ifCNV, maf);
        
        int binNum = this.getBinNumAll();
        int pavcnvNum = this.getNumPAVAndCNV(pavId, cnvId);
        System.out.println("There are " + pavId.length + " pav bins " + String.valueOf(pavId.length/(double)binNum));
        System.out.println("There are " + cnvId.length + " cnv bins " + String.valueOf(cnvId.length/(double)binNum));
        System.out.println("There are " + pavcnvNum + " pav bins " + String.valueOf(pavcnvNum/(double)binNum));
    }
   
    private int getNumPAVAndCNV (String[] pavId, String[] cnvId) {
        TreeSet<String> idSet = new TreeSet();
        for (int i = 0; i < pavId.length; i++) {
            idSet.add(pavId[i]);
        }
        for (int i = 0; i < cnvId.length; i++) {
            idSet.add(cnvId[i]);
        }
        return idSet.size();
    }
    
    public double[] getPAtFDR (double fdr) {
        ArrayList<Double> pPAVList = new ArrayList();
        ArrayList<Double> pCNVList = new ArrayList();
        for (int i = 0; i < pMatrix.length; i++) {
            for (int j = 0; j < pMatrix[i].length; j++) {
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    if (k == refIndex) continue;
                    if (ratioMatrix[i][j][k] > 0) {
                        pCNVList.add(pMatrix[i][j][k]);
                    }
                    else if (ratioMatrix[i][j][k] < 0) {
                        pPAVList.add(pMatrix[i][j][k]);
                    }
                }
            }
        }
        Double[] pPAV = pPAVList.toArray(new Double[pPAVList.size()]);
        Double[] pCNV = pCNVList.toArray(new Double[pCNVList.size()]);
        Arrays.sort(pPAV);
        Arrays.sort(pCNV);
        double[] pCut = new double[2];
        pCut[0] = pPAV[(int)(pPAV.length*fdr)];
        pCut[1] = pCNV[(int)(pCNV.length*fdr)];
        return pCut;
    }
    
    public double[] getRatioAtFDR (double fdr) {
        ArrayList<Double> rList = new ArrayList ();
        for (int i = 0; i < ratioMatrix.length; i++) {
            for (int j = 0; j < ratioMatrix[i].length; j++) {
                for (int k = 0; k < ratioMatrix[i][j].length; k++) {
                    if (k == refIndex) continue;
                    rList.add(ratioMatrix[i][j][k]);
                }
            }
        }
        Double[] rArray = rList.toArray(new Double[rList.size()]);
        Arrays.sort(rArray);
        double[] ratioFDR = new double[2];
        ratioFDR[0] = rArray[(int)(rArray.length*fdr)];
        ratioFDR[1] = rArray[(int)(rArray.length*(1-fdr))];
        return ratioFDR;
    }
    
    private String[] mkHapMapFile (String hapMapFileS, boolean[][][] ifPositive, double maf) {
        ArrayList<String> idList = new ArrayList();
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(hapMapFileS), 65536);
            StringBuilder sb = new StringBuilder ();
            sb.append("rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode");
            for (int i = 0; i < this.getTaxaNum(); i++) sb.append("\t").append(taxaName[i]);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < ifPositive.length; i++) {
                for (int j = 0; j < ifPositive[i].length; j++) {
                    double mf = 0;
                    for (int k = 0; k < ifPositive[i][j].length; k++) {
                        if (ifPositive[i][j][k]) mf+=1;
                    }
                    mf = mf/ifPositive[i][j].length;
                    if (mf > 0.5) mf = 1-mf;
                    if (mf < maf) continue;
                    
                    sb = new StringBuilder();
                    String ID = String.valueOf(i+1)+"-"+String.valueOf(j+1);
                    idList.add(ID);
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
            System.out.println(e.toString());
            System.exit(1);
        }
        String[] id = idList.toArray(new String[idList.size()]);
        return id;
    }
    
    private boolean[][][] discoverCNV (double cutoff) {
        boolean[][][] ifCNV = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifCNV[i] = new boolean[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > cutoff) {
                        ifCNV[i][j][k] = true;
                    }
                    else {
                        ifCNV[i][j][k] = false;
                    }
                }
            }
        }
        return ifCNV;
    }
    
    private boolean[][][] discoverCNV (String cnvCutoffFileS) {
        double[][] cnvCutoff = this.getCutoffValue(cnvCutoffFileS);
        boolean[][][] ifCNV = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifCNV[i] = new boolean[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] > cnvCutoff[i][j]) {
                        ifCNV[i][j][k] = true;
                    }
                    else {
                        ifCNV[i][j][k] = false;
                    }
                }
            }
        }
        return ifCNV;
    }
    
    private boolean[][][] discoverPAV (double cutoff) {
        boolean[][][] ifPAV = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifPAV[i] = new boolean[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    if (ratioMatrix[i][j][k] < cutoff) {
                        ifPAV[i][j][k] = true;
                    }
                    else {
                        ifPAV[i][j][k] = false;
                    }
                }
            }
        }
        return ifPAV;
    }
    
    private boolean[][][] discoverPAV (String pavCutoffFileS) {
        double[][] pavCutoff = this.getCutoffValue(pavCutoffFileS);
        boolean[][][] ifPAV = new boolean[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ifPAV[i] = new boolean[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < chrID.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                for (int k = 0; k < taxaName.length; k++) {
                    System.out.println(pMatrix[i][j][k]+"\t"+ratioMatrix[i][j][k]+"\t"+pavCutoff[i][j]);
                    if (ratioMatrix[i][j][k] < pavCutoff[i][j]) {
                        ifPAV[i][j][k] = true;
                    }
                    else {
                        ifPAV[i][j][k] = false;
                    }
                }
            }
        }
        return ifPAV;
    }
    
    private double[][] getCutoffValue (String cutoffFileS) {
        double[][] cutoff = new double[chrID.length][];
        for (int i = 0; i < cutoff.length; i++) {
            cutoff[i] = new double[binStart[i].length];
        }
        try {
            BufferedReader br = new BufferedReader(new FileReader(cutoffFileS), 65536);
            br.readLine();
            for (int i = 0; i < cutoff.length; i++) {
                for (int j = 0; j < cutoff[i].length; j++) {
                    String[] temp = br.readLine().split("\t");
                    cutoff[i][j] = Double.valueOf(temp[2]);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return cutoff;
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
    
    public RatioMatrixB73Permute getPermutedB73Ratiomatrix (BinByTaxa bbt, int sampleNum, int popSize) {
        int refIndex = 0;
        for (int i = 0; i < bbt.getTaxaNum(); i++) {
            if (bbt.taxaName[i].equals("b73")) {
                refIndex = i;
                break;
            }
        }
        short[][][] nCount = new short[this.getChrNum()][][];
        for (int i = 0; i < nCount.length; i++) {
            nCount[i] = new short[bbt.binStart[i].length][popSize];
        }
        for (int i = 0; i < popSize; i++) {
            short[][] taxaCount = this.getAverCount(bbt, sampleNum, refIndex);
            for (int j = 0; j < taxaCount.length; j++) {
                for (int k = 0; k < taxaCount[j].length; k++) {
                    nCount[j][k][i] = taxaCount[j][k];
                }
            }
        }
        for (int i = 0; i < nCount.length; i++) {
            for (int j = 0; j < nCount[i].length; j++) {
                nCount[i][j][0] = bbt.tagBinCount[i][j][refIndex];
            }
        }
        String[] newTaxa = new String[popSize];
        for (int i = 0; i < newTaxa.length; i++) {
            newTaxa[i] = "c";
        }
        newTaxa[0] = "b73";
        BinByTaxa nbbt = new BinByTaxa(newTaxa, bbt.chrID, bbt.binStart, nCount);
        RatioMatrixB73Permute rm = new RatioMatrixB73Permute (nbbt, "b73");
        return rm;
    }
    
    public void mkZscoreTable (String zscoreTableS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (zscoreTableS), 65536);
            bw.write("Chr\tBinStart\tSE\tMean");
            bw.newLine();
            for (int i = 0; i < this.getChrNum(); i++) {
                for (int j = 0; j < this.getBinNumOnChr(i); j++) {
                    DescriptiveStatistics ds = new DescriptiveStatistics(this.ratioMatrix[i][j]);
                    double se = ds.getStandardDeviation()/Math.sqrt(this.getTaxaNum());
                    double mean = ds.getMean();
                    bw.write(String.valueOf(chrID[i])+"\t"+binStart[i][j]+"\t"+String.valueOf(se)+"\t"+String.valueOf(mean));
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
    
    public void mkPermuteCutoffTable (BinByTaxa bbt, String pavCutoffTable, String cnvCutoffTable, double fdr, int sampleNum, int popSize) {
        int refIndex = 0;
        for (int i = 0; i < bbt.getTaxaNum(); i++) {
            if (bbt.taxaName[i].equals("b73")) {
                refIndex = i;
                break;
            }
        }
        short[][][] nCount = new short[this.getChrNum()][][];
        for (int i = 0; i < nCount.length; i++) {
            nCount[i] = new short[bbt.binStart[i].length][popSize];
        }
        for (int i = 0; i < popSize; i++) {
            short[][] taxaCount = this.getAverCount(bbt, sampleNum, refIndex);
            for (int j = 0; j < taxaCount.length; j++) {
                for (int k = 0; k < taxaCount[j].length; k++) {
                    nCount[j][k][i] = taxaCount[j][k];
                }
            }
        }
        for (int i = 0; i < nCount.length; i++) {
            for (int j = 0; j < nCount[i].length; j++) {
                nCount[i][j][0] = bbt.tagBinCount[i][j][refIndex];
            }
        }
        String[] newTaxa = new String[popSize];
        for (int i = 0; i < newTaxa.length; i++) {
            newTaxa[i] = "c";
        }
        newTaxa[0] = "b73";
        BinByTaxa nbbt = new BinByTaxa(newTaxa, bbt.chrID, bbt.binStart, nCount);
        //nbbt.writeTxtFile("M:/test.txt");
        RatioMatrixB73Permute rm = new RatioMatrixB73Permute (nbbt, "b73");
        double[] pCut = rm.getPAtFDR(fdr);
        System.out.println("PAV ratio cutoff is " + pCut[0]);
        System.out.println("CNV ratio cutoff is " + pCut[1]);
        rm.mkPAVCutoffTable(pavCutoffTable, fdr);
        rm.mkCNVCutoffTable(cnvCutoffTable, fdr);
    }
    
    public short[][] getAverCount (BinByTaxa bbt, int sampleNum, int refIndex) {
        int[] index = new int[sampleNum];
        for (int i = 0; i < sampleNum; i++) {
            index[i] = (int)(Math.random()*bbt.getTaxaNum());
            if (index[i] == refIndex) i--;
        }
        short[][] taxaCount = new short[bbt.getChrNum()][];
        for (int i = 0; i < taxaCount.length; i++) {
            taxaCount[i] = new short[bbt.binStart[i].length];
            for (int j = 0; j < taxaCount[i].length; j++) {
                int sum = 0;
                for (int k = 0; k < index.length; k++) {
                    sum += bbt.tagBinCount[i][j][index[k]];
                }
                taxaCount[i][j] = (short)(sum/sampleNum);
            }
        }
        return taxaCount;
    }
    
    public void mkPAVCutoffTable (String pavCutoffTable, double fdr) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(pavCutoffTable), 65536);
            bw.write("Chr\tBinstart\trCutoff\tpCutoff");
            bw.newLine();
            for (int i = 0; i < chrID.length; i++) {
                for (int j = 0; j < binStart[i].length; j++) {
                    
                    ArrayList<Double> pList = new ArrayList();
                    ArrayList<Double> rList = new ArrayList();
                    for (int k = 0; k < taxaName.length; k++) {
                        if (ratioMatrix[i][j][k] < 0) {
                            pList.add(pMatrix[i][j][k]);
                            rList.add(ratioMatrix[i][j][k]);
                        }
                    }
                    Double[] pArray = pList.toArray(new Double[pList.size()]);
                    Double[] rArray = rList.toArray(new Double[rList.size()]);
                    Arrays.sort(pArray);
                    Arrays.sort(rArray); 
                    double rCutoff;
                    double pCutoff;
                    if (rArray.length == 0) {
                        //System.out.println(pArray.length+"\t"+chrID[i]+"\t"+binStart[i][j]+"\t has no PAV ratio < 0");
                        rCutoff = 0;
                        pCutoff = 0.05;
                    }
                    else {
                        rCutoff = rArray[(int)(rArray.length*fdr)];
                        pCutoff = pArray[(int)(rArray.length*fdr)];
                    }
                    
                    bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(binStart[i][j])+"\t"+String.valueOf(rCutoff)+"\t"+pCutoff);
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
    
    public void mkCNVCutoffTable (String cnvCutoffTable, double fdr) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(cnvCutoffTable), 65536);
            bw.write("Chr\tBinstart\trCutoff\tpCutoff");
            bw.newLine();
            for (int i = 0; i < chrID.length; i++) {
                for (int j = 0; j < binStart[i].length; j++) {
                    
                    ArrayList<Double> pList = new ArrayList();
                    ArrayList<Double> rList = new ArrayList();
                    for (int k = 0; k < taxaName.length; k++) {
                        if (ratioMatrix[i][j][k] > 0) {
                            pList.add(pMatrix[i][j][k]);
                            rList.add(ratioMatrix[i][j][k]);
                        }
                    }
                    Double[] pArray = pList.toArray(new Double[pList.size()]);
                    Double[] rArray = rList.toArray(new Double[rList.size()]);
                    Arrays.sort(pArray);
                    Arrays.sort(rArray);
                    double rCutoff;
                    double pCutoff;
                    if (rArray.length == 0) {
                        //System.out.println(pArray.length+"\t"+chrID[i]+"\t"+binStart[i][j]+"\t has no CNV ratio > 0");
                        rCutoff = 0;
                        pCutoff = 0.05;
                    }
                    else {
                        rCutoff = rArray[(int)(rArray.length*(1-fdr))];
                        pCutoff = pArray[(int)(rArray.length*fdr)];
                    }
                    bw.write(String.valueOf(chrID[i])+"\t"+String.valueOf(binStart[i][j])+"\t"+String.valueOf(rCutoff)+"\t"+pCutoff);
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
  
    public void mkRatioMatrixFile (String outfileS) {
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
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void mkPMatrixFile (String outfileS) {
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
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public int getTaxaNum () {
        return taxaName.length;
    }
    
    
    public void transformMatrix () {
        System.out.println("Transforming ratio matrix");
        for (int i = 0; i <  pMatrix.length; i++) {
            System.out.println("On Chromosome " + String.valueOf(i+1));
            for (int j = 0; j < pMatrix[i].length; j++) {
                double min = 0;
                for (int k = 0; k < pMatrix[i][j].length; k++) {
                    if (Double.isNaN(ratioMatrix[i][j][k]) || ratioMatrix[i][j][k] == 0) {
                        ratioMatrix[i][j][k] = 0;
                    }
                    else {
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
                        ratioMatrix[i][j][k] = Double.NaN;
                    }
                    try {
                        pMatrix[i][j][k] = TestUtils.chiSquareTest(count);
                    }
                    catch (Exception e) {
                        System.out.println(e.toString());
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
