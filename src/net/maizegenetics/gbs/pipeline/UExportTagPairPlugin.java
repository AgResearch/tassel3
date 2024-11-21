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
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.util.ArgsEngine;
import org.apache.log4j.Logger;
import net.maizegenetics.genome.BaseEncoder;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author Jason Wallace
 *
 * This plugin takes the TagPairs file produced by earlier steps in the UNEAK
 * pipeline and converts it to a SAM and/or TOPM alignment file that can be used
 * with the normal (reference-genome based) GBS pipeline. Note that most of the
 * information in the SAM file is bogus and serves only as a placeholder, so _do
 * not_ take anything but the actual sequence at face value.
 *
 * The resulting TOPM file can be plugged directly into the normal
 * (reference-based) GBS pipeline with the TagsToSNPByAlignmentPlugin; make sure
 * to use the -y option since the TBT file generated by UNEAK is stored with
 * bytes, not bits. The SAM file from this plugin must first be converted to a
 * TOPM file with the SAMConverterPlugin and then can be run the same way.
 */
public class UExportTagPairPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTBTToMapInfoPlugin.class);
    private String parentDir = null;
    private UTagPairs tp;
    //Variables for input/output
    String infile = "tagPair.tps", sam = "tagPair.sam", topm_bin = "tagPair.topm.bin", topm_text = "tagPair.topm.txt";
    private boolean toSAM = false, toTOPM_bin = false, toTOPM_text = false;	//Flags for what outputs to export to
    int chrom = 0;	//Default dummy chromosome
    int distance = 1000;	//How far to set the dummy coordinates apart
    //Variables needed for manipulating the SAM-format tags
    private static String[] samTag = {"none", "0", "chr1", "0", "255", "64M", "*", "0", "0", "---", "fff", "*",
        "AS:i:100", "XS:i:0", "XN:i:0", "XM:i:0", "XO:i:0", "XG:i:0", "NM:i:0", "YT:Z:UU"};	//Null SAM-format description of a tag
		/*Note: most of these fields, especially the optional ones at the end, are left unchanged and thus don't mean anything (yet). Currently this
     * is just to preserve compatability; in the future they may be altered to reflect reasonable values. (I think the GBS pipeline just ignores them, though - JasonW.)
     */
    private static byte samNameID = 0, samPosID = 3, samCigarID = 5, samTagS = 9, samTagQ = 10;	//Indexes of needed information to output to the SAM file
    private String samTagQuality = StringUtils.repeat("f", Byte.MAX_VALUE);	//Null quality flag; set to Byte.MAX_VALUE b/c length is currently stored in bytes.
    //Variables needed for exporting to TOPM format; most are just filler
    int topm_maxvariants = 8, topm_chrom = 1;
    byte topm_strand = 1, topm_multimaps = 1;
    byte topm_divergence = Byte.MIN_VALUE, topm_multimap_num = Byte.MIN_VALUE, topm_multimap_pos = Byte.MIN_VALUE, topm_dcoP = Byte.MIN_VALUE, topm_mapP = Byte.MIN_VALUE;
    String variantData = "";

    public UExportTagPairPlugin() {
        super(null, false);
    }

    public UExportTagPairPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        //TODO: I have redundant functionality here. Really need to just create a TOPM object and have it create the output info
        tp = new UTagPairs(infile);

        if (toSAM) {
            convertTagPairToSam();
        }
        if (toTOPM_bin || toTOPM_text) {
            TagsOnPhysicalMap myTopm = makeTopmFromTagPairs(tp);
            if (toTOPM_bin) {
                myTopm.writeBinaryFile(new File(topm_bin));
            }
            if (toTOPM_text) {
                myTopm.writeTextFile(new File(topm_text));
            }
        }
        return null;
    }

    /**
     * Take an existing tagPair file and convert it into a SAM file
     */
    private void convertTagPairToSam() {
        logger.info("Printing SAM file to " + sam);
        File samFile = new File(sam);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(samFile));
            String header = "@SQ\tSN:chr0\tLN:" + ((tp.getTagNum() + 1) * distance / 2) + "\n" + "@PG\tID:bowtie2\tPN:bowtie2\tVN:0.0.0\n";		//Make fake header (2 rows)
            bw.append(header);
            bw.append("@CO\tWARNING: The positions, quality scores, and basically everything herein EXCEPT the sequence are placeholders. DO NOT treat them as real.\n");
            for (int i = 0; i < tp.getTagNum(); i += 2) {
                long pos = (i / 2) * distance + 1;
                bw.append(tagToSamFormat(tp, i, pos));
                bw.append(tagToSamFormat(tp, i + 1, pos));
            }
            bw.close();
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(UExportTagPairPlugin.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private String tagToSamFormat(UTagPairs tp, int i, long pos) {
        String[] mytag = Arrays.copyOf(samTag, samTag.length);
        byte mylength = tp.getTagLength(i);
        mytag[samNameID] = "length=" + mylength + "count=1";
        mytag[samPosID] = Long.toString(pos);
        mytag[samCigarID] = mylength + "M";
        mytag[samTagS] = BaseEncoder.getSequenceFromLong(tp.getTag(i));
        mytag[samTagQ] = samTagQuality.substring(0, mylength - 1);
        return (StringUtils.join(mytag, "\t") + "\n");
    }

    /**
     * Convert the TagPair object to a TOPM object, filling in things with
     * mostly dummy data (mostly happens automatically)
     */
    private TagsOnPhysicalMap makeTopmFromTagPairs(UTagPairs tp) {
        HelperTags tempTags = new HelperTags(tp.getTagNum(), tp.getTag(0).length);
        //Load in all the data for each tag pair; having trouble finding the function to load actual sequences
        for (int i = 0; i < tp.getTagNum(); i++) {
            tempTags.setTag(i, tp.getTag(i), tp.getTagLength(i));
        }
        TagsOnPhysicalMap myTopm = new TagsOnPhysicalMap(tempTags);
        for (int i = 0; i < tp.getTagNum(); i++) {
            int pos = (int) Math.floor(i / 2) * distance + 1;
            myTopm.setChromoPosition(i, topm_chrom, topm_strand, pos, pos + tp.getTagLength(i) - 1);
            myTopm.setDivergence(i, (byte) 0);
            myTopm.setDcoP(i, Byte.MIN_VALUE);
            myTopm.setMapP(i, Byte.MIN_VALUE);
            for (int j = 0; j < myTopm.maxVariants; j++) {
                myTopm.setVariantDef(i, j, Byte.MIN_VALUE);
                myTopm.setVariantPosOff(i, j, Byte.MIN_VALUE);
            }
        }
        return myTopm;
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -d (or --distance) distance to pad tag pairs by (default: 1000)\n"
                + " -i (or --input) tagPair file to convert (required) \n"
                + " -m (or --topm) File to export tag pairs to in binary TOPM format\n"
                + " -s (or --sam) File to output tag pairs to in SAM format\n"
                + " -t (or --topm_text) File to export tag pairs to in TOPM format\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-d", "--distance", true);
            engine.add("-i", "--input", true);
            engine.add("-m", "--topm", true);
            engine.add("-s", "--sam", true);
            engine.add("-t", "--topm_text", true);
            engine.parse(args);
        }
        if (engine.getBoolean("-d")) {
            distance = Integer.parseInt(engine.getString("-d"));
        }
        if (engine.getBoolean("-i")) {
            infile = engine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease specify an input file.\n\n");
        }
        if (engine.getBoolean("-m")) {
            topm_bin = engine.getString("-m");
            toTOPM_bin = true;
        }
        if (engine.getBoolean("-s")) {
            sam = engine.getString("-s");
            toSAM = true;
        }
        if (engine.getBoolean("-t")) {
            topm_text = engine.getString("-t");
            toTOPM_text = true;
        }

        //Check that at least one output format selected
        if (!toSAM && !toTOPM_bin && !toTOPM_text) {
            logger.warn("Warning! Must select either SAM or TOPM format to export to");
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

/*
 * A small helper class that exists solely to pass tag info to TOPM
 */
class HelperTags extends AbstractTags {

    /* Inherited
     protected int tagLengthInLong;  //
     protected long[][] tags;  //Index in rows, actual tag components in columns
     protected byte[] tagLength;  // length of tag (number of bases)  // 1 byte
     */
    public HelperTags(int numtags, int myTagLengthInLong) {
        tagLengthInLong = myTagLengthInLong;
        tags = new long[tagLengthInLong][numtags];
        tagLength = new byte[numtags];
    }

    public void setTag(int index, long[] tagValue, byte myTagLength) {
        tagLength[index] = myTagLength;
        for (int i = 0; i < tagValue.length; i++) {
            tags[i][index] = tagValue[i];
        }
    }
}