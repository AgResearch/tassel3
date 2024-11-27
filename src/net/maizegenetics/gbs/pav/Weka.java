/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

/**
 *
 * @author Fei Lu
 */
public class Weka {
    String[] header;
    int sliceNum = 200;
    File subInputDir;
    File subOutputDir;
    public Weka () {}
    
    public Weka (String inputFileS, String modelFileS, String outputFileS) {
        this.sliceInput(inputFileS);
        this.slicePredict(modelFileS);
        this.mergePredict(outputFileS);
    }
    
    private void mergePredict (String outputFileS) {
        File[] predictFiles = subOutputDir.listFiles();
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outputFileS), 65536);
            bw.write("inst#\tactual\tpredicted\terror");
            bw.newLine();
            int cnt = 1;
            for (int i = 0; i < predictFiles.length; i++) {
                BufferedReader br = new BufferedReader (new FileReader(predictFiles[i]), 65536);
                for (int j = 0; j < 5; j++) br.readLine();
                String temp;
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split("\\s+");
                    if (tem.length < 5) continue;
                    bw.write(String.valueOf(cnt)+"\t"+tem[2]+"\t"+tem[3]+"\t"+tem[4]);
                    bw.newLine();
                    cnt++;
                }
                br.close();
            }
            bw.flush();
            bw.close();
            
        }
        catch (Exception e) {
            System.err.println(e.toString());
    	      e.printStackTrace();
            System.exit(1);
        }
        this.deleteSub(subInputDir);
        this.deleteSub(subOutputDir);
        System.out.println("Prediction files are merged at " + outputFileS);
    }
    private void deleteSub (File dir) {
        File[] files = dir.listFiles();
        for (int i = 0; i < files.length; i++) {
            files[i].delete();
        }
        dir.delete();
    }
    
    private void slicePredict (String modelFileS) {
        File[] inputFiles = subInputDir.listFiles();
        for (int i = 0; i < inputFiles.length; i++) {
            String name = inputFiles[i].getName().replace("arff", "pre.txt");
            name = new File (subOutputDir, name).toString();
            try {
                Runtime run = Runtime.getRuntime();
                String cmd = "cmd /c java weka.classifiers.rules.M5Rules -p 0 -T " + inputFiles[i] + " -l " + modelFileS + " > " +name;
                Process p = run.exec(cmd);
                p.waitFor();
                System.out.println("Prediction is made at " + name);
            }
            catch (Exception e) {
                System.err.println(e.toString());
        	      e.printStackTrace();
                System.exit(1);
            }
        }
    }
    
    private void sliceInput (String inputFileS) {
        File inFile = new File (inputFileS);
        File dir = inFile.getParentFile();
        subInputDir = new File (dir, "tempIn");
        subOutputDir = new File (dir, "tempOut");
        subInputDir.mkdir();
        subOutputDir.mkdir();
        try {
            BufferedReader br = new BufferedReader (new FileReader(inFile), 65536);
            int cnt1 = 1, cnt2 = 0;
            String temp;
            ArrayList<String> headerList = new ArrayList();
            while ((temp = br.readLine()) != "@data") {
                if (temp.equals("@data")) {
                     break;
                }
                else {
                    headerList.add(temp+"\n");
                    cnt1++;
                }      
            }
            headerList.add(temp+"\n");
            header = headerList.toArray(new String[headerList.size()]);
            while ((temp = br.readLine()) != null) cnt2++;
            int cnt = cnt1+cnt2;
            int leftover = cnt2%sliceNum;
            int[] size = new int[sliceNum];
            for (int i = 0; i < sliceNum; i++) size[i] = (cnt2-leftover)/sliceNum;
            for (int i = 0; i < leftover; i++) size[i]++;
            br = new BufferedReader (new FileReader(inFile), 65536);
            for (int i = 0; i < cnt1; i++) br.readLine();
            for (int i = 0; i < sliceNum; i++) {
                int len = 5-String.valueOf(i).length();
                String name = "slice-";
                for (int j = 0; j < len; j++) name = name+"0";
                name = name + String.valueOf(i)+".arff";
                name = new File(subInputDir, name).toString();
                BufferedWriter bw = new BufferedWriter (new FileWriter(name), 65536);
                for (int j = 0; j < header.length; j++) bw.write(header[j]);
                for (int j = 0; j < size[i]; j++) {
                    temp = br.readLine();
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }        
        }
        catch (Exception e) {
            System.out.println(e.toString());
        }
    }
    
}
