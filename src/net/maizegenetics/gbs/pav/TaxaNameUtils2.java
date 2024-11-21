/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author Fei Lu
 */
public class TaxaNameUtils2 {
    String[] taxaName;
    
    public TaxaNameUtils2 (String taxaNameTableS) {
        this.readTaxaName(taxaNameTableS);
    }
    
    public void readTaxaName (String taxaNameTableS) {
        Table t = new Table (taxaNameTableS);
        taxaName = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaName[i] = t.content[i][0];
        }
        Arrays.sort(taxaName);
    }
    
    public void mkPanelTaxaNumTable (String paperKeyFileS, String taxaNumTableS) {//using merged taxaName
        Table t = new Table (paperKeyFileS);
        TreeSet<String> panelSet = new TreeSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            panelSet.add(t.content[i][7]);
        }
        String[] panel = panelSet.toArray(new String[panelSet.size()]);
        TreeSet<String>[] panelSets = new TreeSet[panel.length];
        String[][] panelTaxaName = new String[panel.length][];
        int[] panelCount = new int[panel.length];
        for (int i = 0; i < panelSets.length; i++) {
            panelSets[i] = new TreeSet();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int hit = Arrays.binarySearch(panel, t.content[i][7]);
            panelSets[hit].add(t.content[i][3]);
        }
        for (int i = 0; i < panelSets.length; i++) {
            panelTaxaName[i] = panelSets[i].toArray(new String[panelSets[i].size()]);
            Arrays.sort(panelTaxaName[i]);
        }
        for (int i = 0; i < taxaName.length; i++) {
            for (int j = 0; j < panel.length; j++) {
                int hit = Arrays.binarySearch(panelTaxaName[j], taxaName[i]);
                if (hit < 0) continue;
                panelCount[j]++;
            }
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (taxaNumTableS), 65536);
            bw.write("Panel\tTaxaNum");
            bw.newLine();
            for (int i = 0; i < panel.length; i++) {
                bw.write(panel[i]+"\t"+String.valueOf(panelCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkLower2FullNameFileS (String infileS, String outfileS) {
        Table t = new Table (infileS);
        String[] lowerTaxa = new String[taxaName.length];
        for (int i = 0; i < lowerTaxa.length; i++) {
            lowerTaxa[i] = taxaName[i].split(":")[0].toLowerCase();
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("TaxaName");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                for (int j = 0; j < taxaName.length; j++) {
                    if (t.content[i][0].equals(lowerTaxa[j])) {
                        bw.write(taxaName[j]);
                        bw.newLine();
                        break;
                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkMergeTaxaNameFileS (String outfileS1, String outfileS2) {
        TreeSet<String> nameSet = new TreeSet();
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            if (temp[0].equals("B73(PI550473)")) continue;
            if (temp[0].equals("B73Htrhm")) continue;
            if (temp[0].toLowerCase().equals("blank")) continue;
            if (temp[0].toLowerCase().equals("unknown")) continue;
            nameSet.add(temp[0]);
        }
        try {
            //BufferedWriter bw1 = new BufferedWriter (new FileWriter(outfileS1), 65536);
            //bw1.write("TaxaName");
            //bw1.newLine();
            BufferedWriter bw2 = new BufferedWriter (new FileWriter(outfileS2), 65536);
            bw2.write("TaxaName");
            bw2.newLine();
            String[] name = nameSet.toArray(new String[nameSet.size()]);
            nameSet = new TreeSet ();
            for (int i = 0; i < name.length; i++) {
                //bw1.write(name[i]);
                //bw1.newLine();
                nameSet.add(name[i].toLowerCase());
            }
            //bw1.flush();
            //bw1.close();
            name = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < name.length; i++) {
                bw2.write(name[i]);
                bw2.newLine();
            }
            bw2.flush();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("There is " + nameSet.size() + " taxa in CNV paper");
    }
}
