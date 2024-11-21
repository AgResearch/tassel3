/*
 * Plugin for the UNEAK pipeline to convert the Tag Pair file into a (fake) TOPM file for reintegration with the normal GBS pipeline
 */
package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.gbs.tagdist.UTagPairs;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
import net.maizegenetics.genome.BaseEncoder;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author Jason Wallace This plugin takes the TagPairs file produced by earlier
 * steps in the UNEAK pipeline and converts it to a SAM alignment file that can
 * be used with the normal (reference-genome based) GBS pipeline. Most of the
 * information in the SAM file is bogus and serves only as a placeholder, so _do
 * not_ take anything but the actual sequence at face value.
 *
 * Note: To plug this back into the normal GBS pipeline, run the
 * SAMConverterPlugin to turn the SAM file into a TOPM (Tags On Physical Map)
 * file, and then combine with the TBT (Tags By Taxa) file generated earlier by
 * UNEAK. Note that UNEAK makes a TBT _BYTE_ file, not a BIT file, and so the -y
 * flag will be needed to get the TagsToSNPByAlignmentPlugin to function
 * properly.
 */
public class UTagPairToSAMPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTBTToMapInfoPlugin.class);
    private String parentDir = null;
    //Variables needed for manipulating the SAM-format tags
    private static String[] samTag = {"none", "0", "chr1", "0", "255", "64M", "*", "0", "0", "---", "fff", "*",
        "AS:i:100", "XS:i:0", "XN:i:0", "XM:i:0", "XO:i:0", "XG:i:0", "NM:i:0", "YT:Z:UU"};	//Null SAM-format description of a tag
		/*Note: most of these fields, especially the optional ones at the end, are left unchanged and thus don't mean anything (yet). Currently this
     * is just to preserve compatability; in the future they may be altered to reflect reasonable values. (I think the GBS pipeline just ignores them, though - JasonW.)
     */
    private static byte nameID = 0, posID = 3, cigarID = 5, tagS = 9, tagQ = 10;	//Indexes of needed information to output to the SAM file
    private String tagQuality = StringUtils.repeat("f", Byte.MAX_VALUE);	//Null quality flag; set to Byte.MAX_VALUE b/c length is currently stored in bytes.

    public UTagPairToSAMPlugin() {
        super(null, false);
    }

    public UTagPairToSAMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -w  Working directory to contain subdirectories\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        System.out.println("Outputting TagPairs as SAM file");
        File pd = new File(parentDir);
        String tagPairFileS = new File(pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.tps";
        String samFileS = new File(pd, UCreatWorkingDirPlugin.childDir[4]).getAbsolutePath() + "/tagPair.sam";
        this.createSam(tagPairFileS, samFileS);
        System.out.println("Tag pairs output to " + samFileS);
        return null;
    }

    /**
     * Take an existing tagPair file and convert it into a SAM file
     *
     * @param tagPairFileS The absolute filepath of the tag pair input file as a
     * string
     * @param samFileS The absolute filepath of the SAM output file as a string
     */
    private void createSam(String tagPairFileS, String samFileS) {
        System.out.println("Printing SAM file to " + samFileS);
        UTagPairs tp = new UTagPairs(tagPairFileS);
        File samFile = new File(samFileS);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(samFile));
            String header = "@SQ\tSN:chr0\tLN:" + ((tp.getTagNum() + 1) * 500) + "\n" + "@PG\tID:bowtie2\tPN:bowtie2\tVN:0.0.0\n";		//Make fake header (2 rows)
            bw.append(header);
            for (int i = 0; i < tp.getTagNum(); i += 2) {
                long pos = i * 1000 + 1;
                bw.append(tagToSamFormat(tp, i, pos));
                bw.append(tagToSamFormat(tp, i + 1, pos));
            }
            bw.close();
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(UTagPairToSAMPlugin.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private String tagToSamFormat(UTagPairs tp, int i, long pos) {
        String[] mytag = Arrays.copyOf(samTag, samTag.length);
        byte mylength = tp.getTagLength(i);
        mytag[nameID] = "length=" + mylength + "count=1";
        mytag[posID] = Long.toString(pos);
        mytag[cigarID] = mylength + "M";
        mytag[tagS] = BaseEncoder.getSequenceFromLong(tp.getTag(i));
        mytag[tagQ] = tagQuality.substring(0, mylength - 1);
        return (StringUtils.join(mytag, "\t") + "\n");
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
            engine.parse(args);
        }
        if (engine.getBoolean("-w")) {
            parentDir = engine.getString("-w");
        } else {
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
