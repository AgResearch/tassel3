/*
 * TextToBinaryPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.UTagPairs;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class UTextToBinaryPlugin extends AbstractPlugin {

    private Logger myLogger = Logger.getLogger(UTextToBinaryPlugin.class);

    private ArgsEngine myEngine = null;
    private String myInput;
    private String myOutput;
    private UBinaryToTextPlugin.FILE_TYPES myType;

    public UTextToBinaryPlugin () {
		super(null, false);
	}
    
    public UTextToBinaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        switch (getType()) {
            case TagPairs:
                UTagPairs tps = new UTagPairs ();
                tps.readTxtTagPair(myInput);
                tps.writeTagPair(myOutput);
                break;
        }

        return null;

    }

    @Override
    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myEngine == null) {
            myEngine = new ArgsEngine();
            myEngine.add("-i", "--input-file", true);
            myEngine.add("-o", "--output-file", true);
            myEngine.add("-t", "--file-type", true);
        }

        myEngine.parse(args);

        if (myEngine.getBoolean("-i")) {
            myInput = myEngine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the input file name.");
        }

        if (myEngine.getBoolean("-o")) {
            myOutput = myEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the output file name.");
        }

        if (myEngine.getBoolean("-t")) {
            String temp = myEngine.getString("-t");
            if (temp.equalsIgnoreCase(UBinaryToTextPlugin.FILE_TYPES.TagPairs.toString())) {
                setType(UBinaryToTextPlugin.FILE_TYPES.TagPairs);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the file type.");
        }

    }

    private void printUsage() {
        myLogger.info(
                "\nUsage:\n"
                + "TextToBinaryPlugin <options>\n"
                + " -i  Input File Name\n"
                + " -o  Output File Name\n"
                + " -t  File Type (TOPM, TagCounts, TBTBit, TBTByte)\n");
    }

    public void setInput(String filename) {
        myInput = filename;
    }

    public String getInput() {
        return myInput;
    }

    public void setOutput(String filename) {
        myOutput = filename;
    }

    public String getOutput() {
        return myOutput;
    }

    public void setType(UBinaryToTextPlugin.FILE_TYPES type) {
        myType = type;
    }

    public UBinaryToTextPlugin.FILE_TYPES getType() {
        return myType;
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
