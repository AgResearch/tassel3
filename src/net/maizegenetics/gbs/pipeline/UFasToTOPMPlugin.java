/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import javax.swing.ImageIcon;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;

/**
 *
 * @author Jon Zhang
 */
public class UFasToTOPMPlugin extends AbstractPlugin {
    
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTagCountToTagPairPlugin.class);
    private String parentDir = null;
    
    int numTags = 0;

    public UFasToTOPMPlugin () {
        super(null, false);
    }

    public UFasToTOPMPlugin (Frame parentFrame) {
        super(parentFrame, false);
    }
    
    private void printUsage() {
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -w  Working directory to contain subdirectories\n");
    }
    
    // Count number of lines in .fas file, count tags
    private void countLines(String inFile) {
        try {
            System.out.println("Counting number of tags in " + inFile);
            InputStream is = new BufferedInputStream(new FileInputStream(inFile));
            byte[] c = new byte[1024];
            int count = 0;
            int readChars = 0;
            boolean empty = true;
            while ((readChars = is.read(c)) != -1) {
                empty = false;
                for (int i = 0; i < readChars; i++) {
                    if (c[i] == '\n' || c[i] == '\r') {
                        count++;
                    }
                }
            }
            int lines = (count == 0 && !empty) ? 1 : count;
            numTags = lines/2;
            System.out.println(numTags + " found!");
            is.close();
        } catch(Exception e) {
            System.out.println("Error in reading HapMap.fas file: " + e);
        }
    }
    
    // Reading .fas file and writing to .topm file at some time
    private void fasToTOPM(String inFile, String outFile) {
        try {
            System.out.println("Writing TOPM file " + outFile);
            BufferedReader br = new BufferedReader(new FileReader(inFile));
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile)));
            dos.writeBytes(numTags + "\t2\t0\n");
            for (int i = 0; i < numTags; i++) {
                String[] line1 = br.readLine().split("_");
                String line2 = br.readLine();
                
                int length = Integer.parseInt(line1[2]);
                
                StringBuilder sb = new StringBuilder();
                sb.append(sb);
                sb.append(line2 + "\t");
                sb.append(length + "\t");
                sb.append("1\t");
                sb.append("1\t");
                sb.append("1\t");
                sb.append((int)Math.floor(i/2.0)*100 + "\t");
                sb.append((int)(Math.floor(i/2.0)*100 + length - 1) + "\t");
                sb.append("*\t");
//                for (int j = 0; j < maxVariants; j++) {
//                    sb.append(printWithMissing(variantPosOff[j][row])+"\t");
//                    sb.append(printWithMissing(variantDef[j][row])+"\t");
//                }
                sb.append("*\t");
                sb.append("*\t\n");
                dos.writeBytes(sb.toString());
            }
            System.out.println("TOPM file created!");
            br.close();
            dos.close();
        } catch(Exception e) {
            System.out.println("Error in converting to TOPM: " + e);
        }
    }
    
    @Override
    public DataSet performFunction(DataSet input) {
        File pd = new File (parentDir);
        String fasFile = new File (pd, UCreatWorkingDirPlugin.childDir[7]).getAbsolutePath() + "/HapMap.fas.txt";
        String topmFile = new File (pd, UCreatWorkingDirPlugin.childDir[6]).getAbsolutePath() + "/fas.topm.txt";
        try {
            countLines(fasFile);
            fasToTOPM(fasFile, topmFile);
        } catch(Exception e) {
            System.out.println("Error in file convert: " + e);
        }
        return null;
    }
    
    @Override
    public void setParameters(String args[]) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if(engine == null){
            engine = new ArgsEngine();
            engine.add("-w", "--working-directory", true);
            engine.parse(args);
        }
        if (engine.getBoolean("-w")) {
            parentDir = engine.getString("-w");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("Please specify the working directory.");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
