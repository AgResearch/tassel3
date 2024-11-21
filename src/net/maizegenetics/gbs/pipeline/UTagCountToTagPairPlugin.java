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
public class UTagCountToTagPairPlugin extends AbstractPlugin {

	private ArgsEngine engine = null;
	private Logger logger = Logger.getLogger(UTagCountToTagPairPlugin.class);
	private String parentDir = null;
	private String infile, outfile;
	private boolean useWorkDir = true;
	private double etr = 0.03;

	public UTagCountToTagPairPlugin() {
		super(null, false);
	}

	public UTagCountToTagPairPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}

	private void printUsage() {
		logger.info(
				"\n\nUsage is as follows:\n"
				+ " -e  Error tolerance rate in the network filter. (Default: " + etr + ")\n"
				+ " -i  Input file of merged tag counts\n"
				+ " -o  Output file of tag pairs\n"
				+ " -w  Working directory to contain subdirectories\n"
				+ "(note: must supply EITHER the working directory OR both the input and output files\n)" );
	}

	@Override
	public DataSet performFunction(DataSet input) {
		File pd;
		String mergedTagCountOfAllS, tagPairS;
		if (useWorkDir) {	//if just specified working directory, use defaults (=traditional, maintained for backwards-compatibility)
			pd = new File(parentDir);
			mergedTagCountOfAllS = new File(pd, UCreatWorkingDirPlugin.childDir[3]).getAbsolutePath() + "/mergedAll.cnt";
			tagPairS = new File(pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.tps";
		} else {	//If specific input/output files specified, use them instead.
			mergedTagCountOfAllS = new File(infile).getAbsolutePath();
			tagPairS = new File(outfile).getAbsolutePath();
		}
		UNetworkFilter unf = new UNetworkFilter(mergedTagCountOfAllS, etr, tagPairS);
		return null;
	}

	@Override
	public void setParameters(String[] args) {
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}
		if (engine == null) {
			engine = new ArgsEngine();
			engine.add("-w", "--working-directory", true);
			engine.add("-e", "--error-tolerance-rate", true);
			engine.add("-i", "--input", true);
			engine.add("-o", "--output", true);
			engine.parse(args);
		}
		if (engine.getBoolean("-w")) {
			parentDir = engine.getString("-w");
			useWorkDir = true;
		} else {
			if (engine.getBoolean("-i") && engine.getBoolean("-o")) {
				infile = engine.getString("-i");
				outfile = engine.getString("-o");
				useWorkDir = false;
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify either the working directory or both the input and output files.");
			}
		}
		if (engine.getBoolean("-e")) {
			etr = Double.parseDouble(engine.getString("-e"));
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
