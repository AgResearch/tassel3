/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class SmoothRatioMatrix {
    String[] taxaName;
    int[] chrID;
    int[][] binStart;
    double[][][] ratioMatrix;
    public SmoothRatioMatrix (String genomeInfoFileS, int binSize, RatioMatrixVariableBin rm) {
        this.creatBins(genomeInfoFileS, binSize);
        this.calculate(rm, binSize);
    }
    
    private void calculate (RatioMatrixVariableBin rm, int binSize) {
        taxaName = rm.taxaName;
        ratioMatrix = new double[chrID.length][][];
        for (int i = 0; i < chrID.length; i++) {
            ratioMatrix[i] = new double[binStart[i].length][taxaName.length];
        }
        for (int i = 0; i < binStart.length; i++) {
            for (int j = 0; j < binStart[i].length; j++) {
                ratioMatrix[i][j] = this.getSmoothValue(i, binStart[i][j], binStart[i][j]+binSize-1, rm);
            }
        }
    }
    
    private double[] getSmoothValue (int chrIndex, int binS, int binE, RatioMatrixVariableBin rm) {
        int indexS = Arrays.binarySearch(rm.binStart[chrIndex], binS);
        if (indexS < 0) indexS = -indexS-2;
        int indexE = Arrays.binarySearch(rm.binStart[chrIndex], binE);
        if (indexE < 0) indexE = -indexE-2;
        int[] overlap = new int[indexE-indexS+1];
        System.out.println(binS);
        for (int i = indexS; i < indexE+1; i++) {
            int rmStart = rm.binStart[chrIndex][i];
            int rmEnd = rm.binStart[chrIndex][i]+rm.binSize[chrIndex][i];
            if (rmStart <= binS) {
                if ((rmEnd) <= binE) {
                    overlap[i-indexS] = rmEnd-binS;
                }
                else {
                    overlap[i-indexS] = binE-binS;
                }
            }
            else {
                if ((rmEnd <= binE)) {
                    overlap[i-indexS] = rmEnd-rm.binStart[chrIndex][i];
                }
                else {
                    overlap[i-indexS] = binE- rm.binStart[chrIndex][i];
                }
            }
        }
        double[] portion = new double[overlap.length];
        int sum = 0;
        for (int i = 0; i < overlap.length; i++) sum+=overlap[i];
        for (int i = 0; i < portion.length; i++) portion[i] = (double)overlap[i]/sum;
        double[] value = new double[rm.taxaName.length];
        for (int i = 0; i < value.length; i++) {
            for (int j = indexS; j < indexS+portion.length; j++) {
                value[i] += Math.pow(2, rm.ratioMatrix[chrIndex][j][i]) * portion[j-indexS];
            }
            value[i] = Math.log(value[i])/Math.log(2);
        }
        return value;
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
    
    public void writeRatioFileS (String outfileS) {
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
}
