/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author Fei Lu
 */
public class TaxaNameUtils {
    String[] allTaxaName;
    
    public TaxaNameUtils (String allTaxaNameFileS) {
        this.readAllTaxa(allTaxaNameFileS);
        Arrays.sort(allTaxaName);
    }
    
    public void readAllTaxa (String allTaxaNameFileS) {
        Table t = new Table (allTaxaNameFileS);
        allTaxaName = this.getTaxaName(t);
    }
    
    public void mk282List (String panel282FileS, String outfile282FileS) {
        Table t = new Table (panel282FileS);
        String[] panelTaxaName = this.getTaxaName(t);
        panelTaxaName = this.getLowerName(panelTaxaName);
        panelTaxaName = this.getNonredudentName(panelTaxaName);
        String[] commonName = this.getCommonName(panelTaxaName);
        this.writeTaxaNamelist(commonName, outfile282FileS);
    }
    
    public void mkAmesList (String panelAmesFileS, String outfileAmesFileS) {
        Table t = new Table (panelAmesFileS);
        String[] panelTaxaName = this.getTaxaName(t);
        panelTaxaName = this.getGeneralName(panelTaxaName);
        panelTaxaName = this.getLowerName(panelTaxaName);
        panelTaxaName = this.getNonredudentName(panelTaxaName);
        String[] commonName = this.getCommonName(panelTaxaName);
        this.writeTaxaNamelist(commonName, outfileAmesFileS);
    }
    
    public void mkNAMParentsFullName (String namParentFileS, String namParentFullnameFileS) {
        Table t = new Table (namParentFileS);
        String[] name = new String[t.getRowNumber()];
        for (int i = 0; i < name.length; i++) name[i] = t.content[i][0].toLowerCase();
        Arrays.sort(name);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (namParentFullnameFileS), 65536);
            bw.write("TaxaName");
            bw.newLine();
            for (int i = 0; i < this.allTaxaName.length; i++) {
                String query = allTaxaName[i].split(":")[0].toLowerCase();
                int hit = Arrays.binarySearch(name, query);
                if (hit < 0) continue;
                bw.write(this.allTaxaName[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkB73InitialList (String b73FileS) {
        ArrayList<Integer> indexList = new ArrayList();
        for (int i = 0; i < allTaxaName.length; i++) {
            if (allTaxaName[i].split(":")[0].equals("B73(PI550473)")) indexList.add(i);
            if (allTaxaName[i].split(":")[0].equals("B73")) indexList.add(i);
            //if (allTaxaName[i].split(":")[0].equals("B73Htrhm")) indexList.add(i);
            if (allTaxaName[i].split(":")[0].equals("b73")) indexList.add(i);
        }
        String[] taxaName = new String[indexList.size()];
        for (int i = 0; i < taxaName.length; i++) {
            taxaName[i] = allTaxaName[indexList.get(i)];
        }
        this.writeTaxaNamelist(taxaName, b73FileS);
    }
    
    public void mk282AmesList (String panel282FileS, String panelAmesFileS, String out282AmesFileS) {
        Table t = new Table (panel282FileS);
        String[] panel282TaxaName = this.getTaxaName(t);
        Table t2 = new Table (panelAmesFileS);
        String[] panelAmesTaxaName = this.getTaxaName(t2);
        panelAmesTaxaName = this.getGeneralName(panelAmesTaxaName);
        String[] taxaName = new String[panel282TaxaName.length+panelAmesTaxaName.length];
        System.arraycopy(panel282TaxaName, 0, taxaName, 0, panel282TaxaName.length);
        System.arraycopy(panelAmesTaxaName, 0, taxaName, panel282TaxaName.length, panelAmesTaxaName.length);
        taxaName = this.getNonredudentName(taxaName);
        taxaName = this.getFullNameWithb73(taxaName);
        this.writeTaxaNamelist(taxaName, out282AmesFileS);
    }
    
    public void writeTaxaNamelist (String[] commonName, String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS),  65536);
            bw.write("TaxaName");
            bw.newLine();
            int hit = Arrays.binarySearch(commonName, "b73");
            if (hit < 0) {
                bw.write("b73");
                bw.newLine();
            }
            for (int i = 0; i <  commonName.length; i++) {
                bw.write(commonName[i]);
                bw.newLine();
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
    
    private String[] getTaxaName (Table t) {
        String[] taxaName = new String[t.getRowNumber()];
        for (int i = 0; i < taxaName.length; i++) {
            taxaName[i] = t.content[i][0];
        }
        Arrays.sort(taxaName);
        return taxaName;
    }
    
    private String[] getNonredudentName (String[] taxaName) {
        TreeSet<String> nameSet = new TreeSet();
        for (int i = 0; i < taxaName.length; i++) nameSet.add(taxaName[i]);
        String[] name = nameSet.toArray(new String[nameSet.size()]);
        return name;
    }
    
    private String[] getGeneralName (String[] taxaName) {
        for (int i = 0; i < taxaName.length; i++) {
            String[] temp = taxaName[i].split(":");
            taxaName[i] = temp[0];
        }
        return taxaName;
    }
    
    private String[] getFullNameWithb73 (String[] taxaName) {
        Arrays.sort(taxaName);
        boolean[] iffind = new boolean[taxaName.length];
        ArrayList<String> nameList = new ArrayList();
        for (int i = 0; i < allTaxaName.length; i++) {
            String query = allTaxaName[i].split(":")[0];
            int hit = Arrays.binarySearch(taxaName, query);
            if (hit < 0) continue;
            if (iffind[hit] == true) continue;
            iffind[hit] = true;
            nameList.add(allTaxaName[i]);
        }
        nameList.add("b73");
        String[] fullName = nameList.toArray(new String[nameList.size()]);
        return fullName;
    }
    
    private String[] getCommonName (String[] taxaName) {
        Arrays.sort(taxaName);
        ArrayList<String> nameList = new ArrayList();
        for (int i = 0; i < allTaxaName.length; i++) {
            int hit = Arrays.binarySearch(taxaName, allTaxaName[i]);
            if (hit < 0) continue;
            nameList.add(allTaxaName[i]);
        }
        String[] commonName = nameList.toArray(new String[nameList.size()]);
        return commonName;
    }
    
    private String[] getLowerName (String[] taxaName) {
        for (int i = 0; i < taxaName.length; i++) {
            taxaName[i] = taxaName[i].toLowerCase();
        }
        return taxaName;
    }
}
