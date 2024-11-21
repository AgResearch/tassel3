/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

/**
 * Used to order UNENK SNPs for switchgrass GWAS based on switchgrass genome 1.0
 * @author fl262
 */
public class GoSwitch {

    public GoSwitch() {
        String originalFastaS = "M:/switchgrass_gwas/genome/Panicum_virgatum.main_genome.scaffolds.fasta";
        String statisticsFileS = "M:/switchgrass_gwas/genome/statistics.txt";
        //this.mkStatistics(originalFastaS, statisticsFileS);

        String formatFastaS = "M:/switchgrass_gwas/genome/switch_format.fasta";
        //this.changeFormat(originalFastaS, formatFastaS);
        
        String samFileS = "M:/switchgrass_gwas/alignment/asso_0.05_0.5_0.5_min5_e3.sam";
        String hapMapFileS = "M:/Generator/HapMap/asso_0.05_0.5_0.5_min5_e3/HapMap.hmp.txt";
        String newHapMapFileS = "M:/switchgrass_gwas/hapMap/asso_0.05_0.5_0.5_min5_e3.hapmap.txt";
        this.assignPosition(samFileS, hapMapFileS, newHapMapFileS);
    }

    public static void main(String[] args) {
        new GoSwitch();
    }

    public void assignPosition (String samFileS, String hapMapFileS, String newHapMapFileS) {
        ArrayList<String> nameList = new ArrayList();
        TreeMap<String, String> posMap = new TreeMap();
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(samFileS), 65536);
            String temp;
            while (!(temp = br.readLine()).startsWith("@PG")) {}
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (tem[1].equals("4")) {
                    br.readLine();
                    continue;
                }
                String query = tem[0];
                String strandQ;
                if (tem[1].equals("0")) {
                    strandQ = "+";
                }
                else {
                    strandQ = "-";
                }
                int chrQ = Integer.valueOf(tem[2]);
                int posQ = Integer.valueOf(tem[3]);
                
                temp = br.readLine();
                tem = temp.split("\t");
                if (tem[1].equals("4")) {
                    continue;
                }
                String hit = tem[0];
                String strandH;
                if (tem[1].equals("0")) {
                    strandH = "+";
                }
                else {
                    strandH = "-";
                }
                int chrH = Integer.valueOf(tem[2]);
                int posH = Integer.valueOf(tem[3]);
                if (!strandQ.equals(strandH)) continue;
                if (chrQ != chrH) continue;
                if (posQ != posH) continue;
                String posS = strandH+"_"+String.valueOf(chrH)+"_"+String.valueOf(posH);
                String markerName = hit.split("_")[0];
                nameList.add(markerName);
                posMap.put(markerName, posS);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] markername = nameList.toArray(new String[nameList.size()]);
        Arrays.sort(markername);
        try {
            BufferedReader br = new BufferedReader (new FileReader(hapMapFileS), 65536);
            BufferedWriter bw = new BufferedWriter(new FileWriter(newHapMapFileS), 65536);
            bw.write(br.readLine());
            bw.newLine();
            String temp;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int hit = Arrays.binarySearch(markername, tem[0]);
                if (hit < 0) continue;
                String[] p = posMap.get(tem[0]).split("_");
                tem[2] = p[1];
                tem[3] = p[2];
                tem[4] = p[0];
                for (int i = 0; i < tem.length-1; i++) {
                    bw.write(tem[i]+"\t");
                }
                bw.write(tem[tem.length-1]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void changeFormat(String originalFastaS, String formatFastaS) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(originalFastaS), 65536);
            BufferedWriter bw = new BufferedWriter(new FileWriter(formatFastaS), 65536);
            String temp = null;

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">contig")) {
                    break;
                }
                if (temp.startsWith(">")) {
                    String[] tem = temp.split("_");
                    int chr = Integer.valueOf(tem[0].substring(5, 6));
                    if (tem[0].endsWith("b")) {
                        chr = chr + 9;
                    }
                    temp = ">" + String.valueOf(chr);
                }
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkStatistics(String originalFastaS, String statisticsFileS) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(originalFastaS), 65536);
            BufferedWriter bw = new BufferedWriter(new FileWriter(statisticsFileS), 65536);
            String temp = null;
            int cnt = 0;
            int length = 0;
            String name = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    if (cnt == 0) {
                    } else {
                        bw.write(name + "\t" + String.valueOf(length));
                        bw.newLine();
                    }
                    name = temp;
                    length = 0;
                    cnt++;
                } else {
                    length += temp.length();
                }
            }
            bw.write(name + "\t" + String.valueOf(length));
            bw.newLine();
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
