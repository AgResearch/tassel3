/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class TagAgainstAnchorLoopPlugin extends AbstractPlugin {
    private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTBTToMapInfoPlugin.class);
    private String tbtDirS = null;
    private String anchorDirS = null;
    private String outfileDirS = null;
    private double pThresh = 0.000001;
    private int miniCount = 30;
    private int coreNum = -1;
    
    public TagAgainstAnchorLoopPlugin () {
		super(null, false);
	}

	public TagAgainstAnchorLoopPlugin (Frame parentFrame) {
		super(parentFrame, false);
    }

	private void printUsage (){
        logger.info(
             "\n\nUsage is as follows:\n"
			+ " -t      Directory of TBT files\n"
			+ " -a      Directory of anchor HapMap (Each HapMap file represents one chromosome)\n"
			+ " -o      Directory of output\n"
			+ " -p      P-value threshold. Default: 0.000001\n"
			+ " -m      Minimum taxa count with tag. Default: 30\n"
            + " -c      Core number (1 thread / core). Default: -1 (all cores, 4 thread / core)\n"
        );
    }
    
    public DataSet performFunction (DataSet input) {
		new TagAgainstAnchorLoop(tbtDirS, anchorDirS, outfileDirS, pThresh, miniCount, coreNum);
		return null;
	}
    
    @Override
	public void setParameters(String [] args) {
		if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
		if(engine == null){
			engine = new ArgsEngine();
            engine.add("-t", "--TBT-directory", true);
			engine.add("-a", "--anchor-directory", true);
			engine.add("-o", "--output-directory", true);
			engine.add("-p", "--p-threshold", true);
			engine.add("-m", "--mininum-taxaCount", true);
            engine.add("-c", "--core-number", true);
            engine.parse(args);
        }
		if (engine.getBoolean("-t")) { tbtDirS = engine.getString("-t");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the TBT directory."); }
		if (engine.getBoolean("-a")) { anchorDirS = engine.getString("-a");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the anchor HapMap directory."); }
		if (engine.getBoolean("-o")) { outfileDirS = engine.getString("-o");}
        else{ printUsage(); throw new IllegalArgumentException("Please specify the output directory.");}
        if (engine.getBoolean("-p")) { pThresh = Double.parseDouble(engine.getString("-p"));}
        if (engine.getBoolean("-m")) { miniCount = Integer.parseInt(engine.getString("-m"));}
        if (engine.getBoolean("-c")) { coreNum = Integer.parseInt(engine.getString("-c"));}
    }
	
	public ImageIcon getIcon() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	public String getButtonName() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	public String getToolTipText() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
