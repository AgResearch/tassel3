/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.util.ArgsEngine;
import net.maizegenetics.gbs.util.MutableSimpleAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import org.apache.log4j.Logger;

/**
 * Tools for characterizing and correcting SNPs segregating in bi-parental populations.
 *
 * Error rates are bounded away from zero, but adding 0.5 error to all error rates that
 * that were observed to be zero.
 *
 * @author edbuckler
 */
public class BiParentalErrorCorrectionPlugin extends AbstractPlugin {
    private double[] errorRate;
    private short[][][] popFreq;
    private double[][] popP;
    private double[][] ldByPop;
    private short[] popOfTaxa;
    private double expSegregation=0.5;
    private int minCountForDetectingError=19;
    private double maxErrorRate=0.05;
    private double minMedianLDR2=-1.0;
    private boolean myRemoveUntestedError=true;
    private boolean myRemoveUntestedLD=false;
    private String outHapMap=null;
//    private String outFreqBin=null;  // currently not used
//    private String outErrorRate=null;  // currently not used
    private ArrayList<String> popNames;
    private ArrayList<String> suppliedPopPrefixes=new ArrayList<String>();
    private File pedigreeFile;
    private Alignment a;
    private static ArgsEngine engine = new ArgsEngine();
    int start=1, end=1;
    private String infile;
    private static final Logger myLogger = Logger.getLogger(QseqToTBTPlugin.class);
    Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
    int homoweight=1;  //weight of homozygous genotypes to heterozgyous genotypes
    //1 assumes low coverage and each homozygote is only counted once, while hets get a count for each allele
    //2 assumes better coverage, and homozygotes are counted twice.
    private double minDistortionRatio=2.5;  //this is minimum distortion to be counted as an error
      //e.g. 2.0, means that the allele frequency has to be less the 25% if the expected is 50%
      

    public BiParentalErrorCorrectionPlugin(){
        super(null, false);
    }

    public BiParentalErrorCorrectionPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }


    private void removeErrorsInAlignment(double relRatio) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte mjb=a.getMajorAllele(s);
            byte mnb=a.getMinorAllele(s);
            byte hetb=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(mjb, mnb);
            byte[] pmnb=new byte[popNames.size()];
            byte[] pmjb=new byte[popNames.size()];
            double[] errVSeg=new double[popNames.size()];
            for (int pop = 0; pop < popNames.size(); pop++) {
                int pmj, pmn;  //pop major and minor

                if(popFreq[0][pop][s]<popFreq[1][pop][s]) {  //major and minor reversed in this population
                    pmj=popFreq[1][pop][s]; pmn=popFreq[0][pop][s]; pmnb[pop]=mjb; pmjb[pop]=mnb;
                } else {
                    pmj=popFreq[0][pop][s]; pmn=popFreq[1][pop][s]; pmnb[pop]=mnb; pmjb[pop]=mjb;
                }
                int present=pmj+pmn;
                if(present<1) {
                    errVSeg[pop]=0;}
                else {
                    binomFunc.setNandP(present,errorRate[s]);
                    double errorP = binomFunc.pdf(pmn);
                    binomFunc.setNandP(present,expSegregation);
                    double segP = binomFunc.pdf(pmn);
                    errVSeg[pop]=errorP/segP;
//                    System.out.println(errorRate[s]+"\t"+pmj+"\t"+pmn+"\t"+errorP+"\t"+segP+"\t"+errVSeg[pop]);
                }
            }
            for (int t = 0; t < a.getSequenceCount(); t++) {
                if(popOfTaxa[t]<0) continue;
                if(errVSeg[popOfTaxa[t]]<1) continue;
                byte b=a.getBase(t, s);
                if(b==DataType.UNKNOWN_BYTE) continue;
                if(b==pmnb[popOfTaxa[t]]) msa.setBase(t, s, DataType.UNKNOWN_BYTE);
                if(b==hetb) msa.setBase(t, s, pmjb[popOfTaxa[t]]);
            }
         }
        a=msa;
    }

    private void calcLDByPop(Alignment a) {
        int minCntForLD=10;
        ldByPop=new double[popNames.size()][a.getSiteCount()];
        IdGroup idg=a.getIdGroup();
        for (int pop = 0; pop < popNames.size(); pop++) {
            Arrays.fill(ldByPop[pop], Double.NaN);
            ArrayList<Identifier> keepNames=new ArrayList<Identifier>();
            for (int i = 0; i < popOfTaxa.length; i++) if(popOfTaxa[i]==pop) keepNames.add(idg.getIdentifier(i));
            Identifier[] ids=new Identifier[keepNames.size()];
            keepNames.toArray(ids);
            Alignment pa = new MutableSimpleAlignment(FilterAlignment.getInstance(a, new SimpleIdGroup(ids)));
            int[] segSites=AlignmentFilterByGBSUtils.getLowHetSNPs(pa, false, -2.0, minCntForLD, 0.15, 2);
            Alignment paf = new MutableSimpleAlignment(FilterAlignment.getInstance(pa, segSites));
            int windowSize=paf.getSiteCount()/20;
            LinkageDisequilibrium theLD=new LinkageDisequilibrium(paf,
                    minCntForLD, windowSize, LinkageDisequilibrium.testDesign.SlidingWindow);
            theLD.run();
            for (int i = 0; i < paf.getSiteCount(); i++) {
                ArrayList<Double> obsR2=new ArrayList<Double>();
                for(int j=i-(windowSize/2); j<i+(windowSize/2); j++) {
                    if(Double.isNaN(theLD.getRSqr(i, j))) continue;
                    if(Math.abs(paf.getPositionInLocus(i)-paf.getPositionInLocus(j))<100000) continue;
                    obsR2.add(theLD.getRSqr(i, j));
                }
                Collections.sort(obsR2);
                if(obsR2.size()>4) {
                    ldByPop[pop][a.getSiteOfPhysicalPosition(paf.getPositionInLocus(i), null)]=obsR2.get(obsR2.size()/2);
                }
            }
  //          System.out.println("POP:"+pop+Arrays.toString(ldByPop[pop]));
        }
    }

    private void removeHighErrorSites(double maxErrorRate, boolean removeUntested) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        int sitesWithHighError=0, untestedSites=0;
        for (int s = 0; s < a.getSiteCount(); s++) {
            if(Double.isNaN(errorRate[s])) {
                if(removeUntested) {
                    msa.clearSite(s);
                    untestedSites++;
                }
            }
            else if(errorRate[s] > maxErrorRate) {
                msa.clearSite(s);
                sitesWithHighError++;
            }
         }
        System.out.printf("Initial Sites %d  ErrorRate %g  HighErrorSites %d UntestedSites %d %n",a.getSiteCount(),
                maxErrorRate, sitesWithHighError, untestedSites);
        msa.sortSiteByPhysicalPosition();
        a=msa;
        System.out.printf("Final Sites %d  ErrorRate %g  HighErrorSites %d %n",a.getSiteCount(), maxErrorRate, sitesWithHighError);
    }

     private void removeLowLDSites(boolean removeUntested) {
        MutableSimpleAlignment msa=new MutableSimpleAlignment(a);
        int sitesWithLowLD=0, untestedSites=0;
        for (int s = 0; s < a.getSiteCount(); s++) {
            ArrayList<Double> obsR2=new ArrayList<Double>();
            for (int pop = 0; pop < popNames.size(); pop++) {
                if(!Double.isNaN(ldByPop[pop][s])) obsR2.add(ldByPop[pop][s]);
            }
            Collections.sort(obsR2);
            double medianR2=Double.NaN;
            if(obsR2.size()>0) medianR2=obsR2.get(obsR2.size()/2);
            if(Double.isNaN(medianR2)) {
                untestedSites++;
                if(removeUntested)  msa.clearSite(s);
            }
            else if(medianR2 < minMedianLDR2) {
                msa.clearSite(s);
                sitesWithLowLD++;
            }
         }
        System.out.printf("Initial Sites %d  minMedianLDR2 %g  LowLDSites %d UntestedSites %d %n",a.getSiteCount(),
                minMedianLDR2, sitesWithLowLD, untestedSites);
        msa.sortSiteByPhysicalPosition();
        a=msa;
        System.out.printf("Final Sites %d  ErrorRate %g  LowLDSites %d %n",a.getSiteCount(), maxErrorRate, sitesWithLowLD);
    }

    private void classifyTaxaToPops(ArrayList<String> suppliedPopPrefixes) {
        popNames=new ArrayList<String>();
        popOfTaxa=new short[a.getSequenceCount()];
        Arrays.fill(popOfTaxa, (short)-1);
        for (int i= 0; i < a.getSequenceCount(); i++) {
            String taxonNameInAlignment=a.getTaxaName(i);
            for (String suppliedPopPrefix : suppliedPopPrefixes) {
                if(taxonNameInAlignment.contains(suppliedPopPrefix)){
                    if(!popNames.contains(suppliedPopPrefix)){
                        popNames.add(suppliedPopPrefix);
                    }
                    popOfTaxa[i]=(short)popNames.indexOf(suppliedPopPrefix);
                }
            }
        }
    }

    private void classifyTaxaToPops(File pedigreeFile) {
        String inputLine = "Nothing has been read from the pedigree input file yet";
        popNames=new ArrayList<String>();
        HashMap<String,Short> popIndexOfTaxaName = new HashMap<String,Short>();
        int famCol = 0, nameCol = 1, contrib1Col = 4, contrib2Col = 5;
        int nTaxa = 0, nTaxaInFamilies = 0;
	try {
	    BufferedReader br = new BufferedReader(new FileReader(pedigreeFile), 65536);
            inputLine = br.readLine();  // skip header line
            while ((inputLine = br.readLine()) != null) {
                ++nTaxa;
                String[] cells = inputLine.split("\t");
                if (cells[famCol].equals("NA") ) continue;
                if (cells[contrib1Col].equals("NA")) continue;
                if (cells[contrib2Col].equals("NA")) continue;
                if (Double.parseDouble(cells[contrib1Col]) != 0.5 ) continue;
                if (Double.parseDouble(cells[contrib2Col]) != 0.5 ) continue;
                if (!popNames.contains(cells[famCol])) popNames.add(cells[famCol]);
                popIndexOfTaxaName.put(cells[nameCol], (short)popNames.indexOf(cells[famCol]));
                ++nTaxaInFamilies;
            }
	} catch (Exception e) {
	    myLogger.error("Catch in reading pedigree file e=" + e);
            e.printStackTrace();
            System.err.println(inputLine);
            System.exit(1);
	}
        myLogger.info(nTaxa + " taxa read from the pedigree file, "+nTaxaInFamilies+" of which belong to a biparental family with expected segregation of 1:1");
        popOfTaxa=new short[a.getSequenceCount()];
        Arrays.fill(popOfTaxa, (short)-1);
        for (int i= 0; i < a.getSequenceCount(); i++) {
            String fullTaxonNameInAlignment=a.getFullTaxaName(i);
            if (popIndexOfTaxaName.containsKey(fullTaxonNameInAlignment)) {
                popOfTaxa[i]=popIndexOfTaxaName.get(fullTaxonNameInAlignment);
            }
        }
    }

    private void classifyTaxaToPops(String popMask) {
        popNames=new ArrayList<String>();
        popOfTaxa=new short[a.getSequenceCount()];
        Arrays.fill(popOfTaxa, (short)-1);
        Pattern pat=Pattern.compile(popMask);
        for (int i=0; i<a.getSequenceCount(); i++) {
            Matcher m=pat.matcher(a.getTaxaName(i));
            if(m.find()) {
                String tn=a.getTaxaName(i).substring(m.start(), m.end());
                if(!popNames.contains(tn)) popNames.add(tn);
                popOfTaxa[i]=(short)popNames.indexOf(tn);
            }
        }
    }

    /**
     * Creates and reports minor allele frequency by population
     * @param a input alignment
     * @param popMask in REGEX format
     */
    private void calcSNPsFreqByPop(Alignment a, double pDevFromExp) {
        popFreq=new short[2][popNames.size()][a.getSiteCount()];
        popP=new double[popNames.size()][a.getSiteCount()];
        errorRate=new double[a.getSiteCount()];
        int[][] errorCorrCnt=new int[2][a.getSiteCount()];
        for (int s = 0; s < a.getSiteCount(); s++) {
            byte mjb=a.getMajorAllele(s);
            byte mnb=a.getMinorAllele(s);
            byte hetb=IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(mjb, mnb);
//            if((mjb==DataType.UNKNOWN_BYTE)||(mnb==DataType.UNKNOWN_BYTE)) continue;
            for (int t = 0; t < a.getSequenceCount(); t++) {
                if(popOfTaxa[t]<0) continue;
                //should this include homo versus heterozygote
                byte b=a.getBase(t, s);
                if(b==DataType.UNKNOWN_BYTE) continue;
                if(b==mjb) popFreq[0][popOfTaxa[t]][s]+=homoweight;
                if(b==mnb) popFreq[1][popOfTaxa[t]][s]+=homoweight;
                if(b==hetb) {popFreq[0][popOfTaxa[t]][s]++; popFreq[1][popOfTaxa[t]][s]++;}
            }
            for (int pop = 0; pop < popNames.size(); pop++) {
                int pmj, pmn;  //pop major and minor
                popP[pop][s]=Double.NaN;
                if(popFreq[0][pop][s]<popFreq[1][pop][s]) {
                    pmj=popFreq[1][pop][s]; pmn=popFreq[0][pop][s];
                } else {
                    pmj=popFreq[0][pop][s]; pmn=popFreq[1][pop][s];
                }
                int present=pmj+pmn;
                if(present<minCountForDetectingError) continue;
                binomFunc.setNandP(present,expSegregation);
                popP[pop][s] = binomFunc.cdf(pmn);
//                System.out.println(pmj+"\t"+pmn+"\t"+popP[pop][s]);
                double minFreq=(double)pmn/(double)present;
                if((popP[pop][s]<pDevFromExp)&&(minFreq<(expSegregation/minDistortionRatio))) {
//                if((popP[pop][s]<pDevFromExp)) {
                    errorCorrCnt[0][s]+=pmn;
                    errorCorrCnt[1][s]+=present;
                }
            }
         }
        for (int s = 0; s < a.getSiteCount(); s++) {
            double error=errorCorrCnt[0][s];
            double total=errorCorrCnt[1][s];
            if(error<1) error=0.5;
            if(total<1) {errorRate[s]=Double.NaN;}  //if the SNP was never present in biparental pops then the error is undefined
            else {errorRate[s]=error/total;}
        }
    }

    public void reportSNPFrequencyByFamily() {
        //TODO make it saveable to file or stdout
        StringBuilder sb=new StringBuilder("Site\tErrorRate\t");
        for (int pop = 0; pop < popNames.size(); pop++) {
                sb.append(popNames.get(pop)); sb.append("\t");
            }
        System.out.println(sb.toString());
        DecimalFormat df=new DecimalFormat("0.00");
        for (int s = 0; s < a.getSiteCount(); s++) {
            sb=new StringBuilder(a.getSNPID(s)+"\t"+errorRate[s]+"\t");
            for (int pop = 0; pop < popNames.size(); pop++) {
                if(Double.isNaN(popP[pop][s])) {sb.append("NaN");}
                else {sb.append(df.format((double)popFreq[1][pop][s]/(double)(popFreq[0][pop][s]+popFreq[1][pop][s])));}
                sb.append("\t");
            }
            System.out.println(sb.toString());
        }
        
    }

    public void reportDistOfFrequency() {
        //TODO make it saveable to file or stdout
        int[] bins=new int[101];
        for (int pop = 0; pop < popNames.size(); pop++) {
            for (int s = 0; s < a.getSiteCount(); s++) {
                if(Double.isNaN(popP[pop][s])) continue;
                double minFreq=(double)popFreq[1][pop][s]/(double)(popFreq[0][pop][s]+popFreq[1][pop][s]);
                bins[(int)Math.round(100*minFreq)]++;
            }
        }
        System.out.println("BinFreqX100\tNumOfSites");
        for (int i = 0; i < bins.length; i++) {
            System.out.printf("%d\t%d%n",i,bins[i]);
        }
    }
    
    public void reportPercentilesOfErrorRates() {
        //TODO make it saveable to file or stdout
        double[] le=Arrays.copyOf(errorRate, errorRate.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i+=(le.length/20)) {
            System.out.printf("%.2g\t%.3g%n",((double)i/(double)le.length),le[i]);
        }
    }

    @Override
    public void setParameters(String[] args){
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-oE", "--outErrorFile", true);
        engine.add("-oB", "--outBinDistFile", true);
        engine.add("-popM", "--popMask", true);
        engine.add("-popF", "--popFile", true);
        engine.add("-pedF", "--pedigreeFile", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-mxE", "--maxError", true);
        engine.add("-mnD", "--minDistortion", true);
        engine.add("-mnPLD", "--minPopLD", true);
        engine.add("-kpUT", "--keepUntested", false);
        engine.parse(args);
        if(engine.getBoolean("-sC")) {
            start=Integer.parseInt(engine.getString("-sC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide a start chromosome.\n");
        }

        if(engine.getString("-popF") != null){
            if(engine.getString("-popM") != null || engine.getString("-pedF") != null) {
                printUsage();
                throw new IllegalArgumentException("Options \"-popF\", \"-popM\", and \"-pedF\"  are mutually exclusive.\n");
            }
            try {
                String currLine;
                BufferedReader br = new BufferedReader(new FileReader(engine.getString("-popF")));
                while((currLine=br.readLine()) != null){
                    suppliedPopPrefixes.add(currLine.trim());
                }
            } catch (Exception e) {
                printUsage();
                throw new IllegalArgumentException("Couldn't open file containing taxon names (option -popF).");
            }
        }
        if(engine.getString("-pedF") != null){
            if(engine.getString("-popM") != null || engine.getString("-popF") != null) {
                printUsage();
                throw new IllegalArgumentException("Options \"-popF\", \"-popM\", and \"-pedF\"  are mutually exclusive.\n");
            }
            pedigreeFile = new File(engine.getString("-pedF"));
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-pedF option: "+engine.getString("-pedF")+")");
            }
        }
        if(engine.getString("-popM") != null){
            if(engine.getString("-popF") != null || engine.getString("-pedF") != null){
                printUsage();
               throw new IllegalArgumentException("Options \"-popF\", \"-popM\", and \"-pedF\"  are mutually exclusive.\n");
            }
        }
        if(engine.getBoolean("-eC")) {
            end=Integer.parseInt(engine.getString("-eC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide an end chromosome.\n");
        }
        infile = engine.getString("-hmp");
        performFunction(null);
    }

//    public static void main(String[] args) {
//            args = new String[]{
//            "-hmp","/media/Data/SAP_070511/SAP_070511.c+.hmp.txt",
//            "-o","/media/Data/SAP_070511/SAP_070511.c+.bpec.hmp.txt",
//            "-oB", "/Users/edbuckler/SolexaAnal/GBS/test/errorBin.txt",
//            "-oE", "/Users/edbuckler/SolexaAnal/GBS/test/errorBySNP.txt",
//            "-popM",".*",
//            "-sC","1","-eC","4",
//            "-mxE","0.01",
//            "-mnD","2.0",
//        };
//
//
//        BiParentalErrorCorrectionPlugin testClass = new BiParentalErrorCorrectionPlugin();
//        testClass.setParameters(args);
//        testClass.performFunction(null);
//    }

//    public void setOutErrorRate(String outErrorRate) {
//        this.outErrorRate = outErrorRate;
//    }

    public void setOutHapMap(String outHapMap) {
        this.outHapMap = outHapMap;
    }

//    public void setOutFreqBin(String outFreqBin) {
//        this.outFreqBin = outFreqBin;
//    }

    public void setMaxErrorRate(double maxErrorRate) {
        this.maxErrorRate = maxErrorRate;
    }

    public double getMinDistortionRatio() {
        return minDistortionRatio;
    }

    public void setMinDistortionRatio(double minDistortionRatio) {
        this.minDistortionRatio = minDistortionRatio;
    }

    public void setMinMedianLDR2(double minMedianLDR2) {
        this.minMedianLDR2 = minMedianLDR2;
    }

        private void printUsage() {
        myLogger.info(
            "\n\n\nAvailable options for the BiParentalErrorCorrectionPlugin are as follows:\n"
            + "-hmp   Input HapMap file\n"
            + "-o     Output HapMap file\n"
//            + "-oE    Output error file\n"
//            + "-oB    Output binary distribution file\n"
            + "One of the following:"
            + "\t  -popM  Population mask: a regular expression specifying the prefixes of each biparental population.  Defaults to all taxa.\n"
            + "\t  -popF  Population file: the name of a file containing the prefixes of each biparental population, one population prefix per line.\n"
            + "\t  -pedF  Pedigree file: lists the population (family) names, the full names, parents, parental contributions\n"
            + "             and the expected F for each taxon in the input HapMap file. Taxa that are not part of a biparental family have \"NA\" as\n"
            + "             their family name.\n"
            + "-sC    Start chromosome\n"
            + "-eC    End chromosome\n"
            +"-mxE    Maximum error\n"
            +"-mnD    Minimum distortion\n"
            +"-mnPLD  Minimum median population LD (R^2)\n"
            +"-kpUT   Keep untested SNPs for error (default remove)\n\n\n"
        );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        while (start <= end) {
            int chr = start;
            ArrayList<Datum> dList=new ArrayList<Datum>();
            String currFile =infile.replace("+", ""+chr);
            System.out.println("Reading: "+currFile);
            Alignment align;
            try {
                align=ImportUtils.readFromHapmap(currFile);
            } catch(Exception e) {
                myLogger.info("Could not read input hapmap file for chr"+chr+":\n\t"+currFile+"\n\tSkipping...");
                continue;
            }
            System.out.println("Finished Reading: "+currFile);
            dList.add(new Datum(currFile,align,"alignment"));
            if(engine.getString("-popM") != null) dList.add(new Datum("PopMask",engine.getString("-popM"),"popmask"));
            if(engine.getBoolean("-o")) setOutHapMap(engine.getString("-o").replace("+", ""+chr));
//            if(engine.getBoolean("-oE")) setOutErrorRate(engine.getString("-oE").replace("+", ""+chr));
//            if(engine.getBoolean("-oB")) setOutFreqBin(engine.getString("-oB").replace("+", ""+chr));
            if(engine.getBoolean("-mxE")) setMaxErrorRate(Double.parseDouble(engine.getString("-mxE")));
            if(engine.getBoolean("-mnD")) setMinDistortionRatio(Double.parseDouble(engine.getString("-mnD")));
            if(engine.getBoolean("-mnPLD")) setMinMedianLDR2(Double.parseDouble(engine.getString("-mnPLD")));
            if(engine.getBoolean("-kpUT")) {this.myRemoveUntestedLD=false; this.myRemoveUntestedError=false;}
            filter(new DataSet(dList,null));
            start++;
        }
        return null;
    }

    private void filter(DataSet input){
        List<Datum> b=input.getDataOfType(Alignment.class);
        Datum c = b.get(0);
        Object d = c.getData();
        a=(Alignment)d;
//        a=(Alignment)input.getDataOfType(Alignment.class).get(0).getData();
        a=new MutableSimpleAlignment(a);
        double realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,false,false);
        double randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,true,false);
        System.out.println("Ratio of RandomToReal:"+randomDist/realDist);

        if(engine.getString("-popF") != null){
            System.out.println("Reading population prefixes from supplied file (-popF option).");
            classifyTaxaToPops(suppliedPopPrefixes);
        } else if (engine.getString("-pedF") != null) {
            System.out.println("Reading population (pedigree) of each individual from supplied pedigree file (-pedF option).");
            classifyTaxaToPops(pedigreeFile);
        } else {
            System.out.println("Reading population prefixes from regex.");
            String pm=(String)input.getDataOfType(String.class).get(0).getData();
            classifyTaxaToPops(pm);
        }

        for (int i = 0; i < suppliedPopPrefixes.size(); i++) {
            System.out.println(suppliedPopPrefixes.get(i));
        }
        for (int i = 0; i < popNames.size(); i++) {
            System.out.println(popNames.get(i));
        }
//        for (int i= 0; i < popOfTaxa.length; i++) {
//            System.out.println(popOfTaxa[i]);
//        }
        if(this.minMedianLDR2>0) {
            calcLDByPop(a);
            removeLowLDSites(myRemoveUntestedLD);
        }
        calcSNPsFreqByPop(a,0.001);
//        reportSNPFrequencyByFamily();
        reportDistOfFrequency();
        reportPercentilesOfErrorRates();
        removeErrorsInAlignment(1.0);
        removeHighErrorSites(maxErrorRate,myRemoveUntestedError);
        calcSNPsFreqByPop(a,0.001);
//        reportSNPFrequencyByFamily();
        reportDistOfFrequency();
        reportPercentilesOfErrorRates();
        realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,false,true);
        randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true,true,false);
        System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
        if(outHapMap!=null) ExportUtils.writeToHapmap(a, false, this.outHapMap, '\t');
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "FilterErrorForBiparental";
    }

    @Override
    public String getToolTipText() {
        return "Filters and esitmates error rates from biparental populations";
    }
}
