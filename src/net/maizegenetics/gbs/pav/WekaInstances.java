/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pav;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

/**
 *
 * @author Fei Lu
 */
public class WekaInstances {
    String[] attNames;
    double[][] value;
    
    public WekaInstances () {}
    
    public WekaInstances (String[] attNames, double[][] value) {
        this.attNames = attNames;
        this.value = value;
    }
    
    public byte[] getCluster (String infileS) {
        byte[] cluster = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
            String temp;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.length() == 0) continue;
                if (temp.startsWith("@")) continue;
                cnt++;
            }
            cluster = new byte[cnt];
            cnt = 0;
            br = new BufferedReader (new FileReader(infileS), 65536);
            while ((temp = br.readLine()) != null) {
                if (temp.length() == 0) continue;
                if (temp.startsWith("@")) continue;
                String[] tem = temp.split(",");
                cluster[cnt] = Byte.valueOf(tem[tem.length-1].replace("cluster", ""));
                cnt++;
            }
            br.close();
        }
        catch (Exception e) {
            System.err.println();
            e.printStackTrace();
            System.exit(1);
        }
        return cluster;
    }
    
    public void writeNumericARFF (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("@relation 'Output from Java'");
            bw.newLine();bw.newLine();
            for (int i = 0; i < attNames.length; i++) {
                 bw.write("@attribute "+attNames[i]+" numeric");
                 bw.newLine();
            }
            bw.newLine();
            bw.write("@data");
            bw.newLine();
            for (int i = 0; i < value[0].length; i++) {
                for (int j = 0; j < value.length-1; j++) {
                    bw.write(String.valueOf(value[j][i])+",");
                }
                bw.write(String.valueOf(value[value.length-1][i]));
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
}
