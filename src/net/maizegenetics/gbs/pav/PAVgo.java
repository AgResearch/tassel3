/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/

package net.maizegenetics.gbs.pav;

import cern.jet.stat.Probability;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesByRowPlugin;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesPlugin;
import net.maizegenetics.gbs.pipeline.TagAgainstAnchor;
import net.maizegenetics.gbs.pipeline.TagAgainstAnchorLongTimePosBlock;
import net.maizegenetics.gbs.pipeline.TagCallerAgainstAnchorMT;
import net.maizegenetics.gbs.pipeline.UFastqToTagCountPlugin;
import net.maizegenetics.gbs.pipeline.UNetworkFilter;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.UTagCountMutable;
import net.maizegenetics.genome.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class PAVgo {
    String parentDir;
    String[] childDir = {"tbtByLane", "mergedTBT", "tagCount", "sliceTBT", "anchor", "mappingResult", "tagAlignment", "pedigree", "PhyGenMapping", "cnv", "train",
        "figurefiles", "prediction", "tags", "translocation", "cnv2"};
    File[] fileOfData;
    
    public PAVgo (String workingDirS) {
        this.parentDir = workingDirS;
        this.creatDir(workingDirS);
        //this.test();
        //this.assoPipe();
        //this.jointLinkagePipe();
        this.PAVPipe();
        //this.cnvPipe();
        //this.cnvPipe2();
        //this.translocationPipe();
        //this.trainPipe();
        //this.figurePipe();
        //this.otherPipe();
    }
    
    public void otherPipe () {
        //this.judyPAV1();
        //this.judyPAV2();
        //this.tiffanyPAV();
        //this.nilPipe();
    }
    
    public void nilPipe () {
        //this.parseFaseq();
        this.judithPAV();
    }
    
    public void judithPAV () {
        String tagCountDirS = "M:/collaboration/judith/tagCounts";
        //String inTagFileS = "M:/pav/judy/chr1_175Mb_190Mb_tag.txt";
        //String outTagFileS = "M:/collaboration/judith/chr1_175Mb_190Mb_tag_NIL.txt";
        
        String inTagFileS = "M:/pav/judy/ht_region_tag.txt";
        String outTagFileS = "M:/collaboration/judith/ht_region_tag_NIL.txt";
        
        PAVUtils util = new PAVUtils();
        util.mkJudyFile3(tagCountDirS, inTagFileS, outTagFileS);
    }
    
    public void parseFaseq () {
        String workingDirS = "M:/collaboration/judith/";
        String arguments = "-w " + workingDirS + " -e ApekI";
        String[] args = arguments.split(" ");
        UFastqToTagCountPlugin uftc = new UFastqToTagCountPlugin();
        uftc.setParameters(args);
        uftc.performFunction(null);
    }
    
    public void tiffanyPAV () {
        String outputTBTFileS = "M:/pav/judy/tbt_NAMParents.byte";
        String taxaNameFileS = "M:/pav/cnv2/taxaName/NAMParents_full.txt";
        String inTagFileS = "M:/pav/judy/chr1_175Mb_190Mb_tag.txt";
        String outTagFileS = "M:/pav/judy/chr1_175Mb_190Mb_tag_NAMParents.txt";
        PAVUtils util = new PAVUtils();
        util.mkJudyFile2(outputTBTFileS, taxaNameFileS, inTagFileS, outTagFileS);
    }
    
    public void judyPAV2 () {
        //Step 1
        String inputTBTFileS = "J:/tbt_20.byte";
        String taxaNameFileS = "M:/pav/cnv2/taxaName/NAMParents_full.txt";
        String outputTBTFileS = "M:/pav/judy/tbt_NAMParents.byte";
        //PAVUtils util = new PAVUtils();
        //util.mkSubTBTByTaxa(inputTBTFileS, taxaNameFileS, outputTBTFileS);
        
        //Step 2
        String inTagFileS = "M:/pav/judy/ht_region_tag.txt";
        String outTagFileS = "M:/pav/judy/ht_region_tag_NAMParents.txt";
        PAVUtils util = new PAVUtils();
        util.mkJudyFile2(outputTBTFileS, taxaNameFileS, inTagFileS, outTagFileS);
    }
    
    public void judyPAV1 () {
        //Step 1
        String inputTBTFileS = "J:/tbt_20.byte";
        String taxaNameFileS = "M:/pav/cnv2/taxaName/TX303.txt";
        String outputTBTFileS = "M:/pav/judy/tbt_TX303.byte";
        //PAVUtils util = new PAVUtils();
        //util.mkSubTBTByTaxa(inputTBTFileS, taxaNameFileS, outputTBTFileS);
        
        //Step 2
        String inTagFileS = "M:/pav/judy/ht_region_tag.txt";
        String outTagFileS = "M:/pav/judy/ht_region_tag_TX303.txt";
        PAVUtils util = new PAVUtils();
        util.mkJudyFile1(outputTBTFileS, inTagFileS, outTagFileS);
    }
    
    public void test () {
        String fileS = "M:/pav/tagAlignment/26Mtag-k2.pmap";
        TagPMap pmap = new TagPMap (fileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter("M:/a.txt"));
            for (int i = 0; i < pmap.getTagCount(); i++) {
                long[] t = pmap.getTag(i);
                bw.write(String.valueOf(i)+"\t"+BaseEncoder.getSequenceFromLong(t));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void translocationPipe () {
        //this.creatUniqueJointMappingFile();
        //this.creatGeneUniqueJointMappingFile();
        //this.creatDiJointMappingFile();
        
    }
    
    public void cnvPipe2 () {//change the binsize, each bin has 100 read count average
        //this.creatBinSizeNAMParents();
        //this.creatVariableBinByTaxaNAM(); //special for the NAM parents (high coverage)
        //this.addRemoveb73VairableBBTNAM();
        //this.mergeTaxaVariableBBTNAM();
        //this.callNAMGenotypeVariableBin();
        
        
        //this.creatBinSize();
        //this.creatVariableBinByTaxa();
        this.selectTaxaPanel();
        //this.taxaStatistics();
        //this.addRemoveb73VairableBBT();
        //this.mergeTaxaVariableBBT();
        //this.selectVariableBBT();
        //this.callPaperGenotypeVariableBin();
        //this.callAmes282GenotypeVariableBin();
        //this.call282GenotypeVariableBin();
    }
    
    public void cnvPipe () {
        //this.creatBinByTaxa();
        
        //this.mergedTaxaInBBT(); //or
        //this.addB73TaxaInBBT(); //for GWAS
        
        
        
        //this.mergeBinsInBBT();
        //this.creatTaxaNameList();
        //this.creatTaxaCoverage();
        //this.selectTaxaInBBT();
        
        
        
        //*************************
        //***option fail
        //this.permutationCutoff();
        //this.creatRatioMatrixByPDis();
        //this.creatGenotypeByPDis();
        //this.creatGenotypeByRDis();
        //*************************
        
        //this.creatRatioMatrixByZscore();//depleted
        
        //this.creatRatioMatrix();
        //this.callPAVGenotypes();
        //this.mkPAVgenotypeFrequency();
        
        
        //this.callCNVgenotypes(); //not deprecated yet, non-mendilian genotype, no LD
        //this.getCNVTagCount();//Deprecated
        //this.getTBTOfCNVTag();//Deprecated
        //this.assignPosToCnvTBTAndMergeTaxa();//Deprecated
        //this.getTaxaTagCountInBin();//Deprecated
    }
    
    public void PAVPipe () {
        //this.mkPredictionValueFile();
        
        //this.mergeHighResPGMap();
        this.creatPAVFPGMap();
        
        //this.sortFPGPaparFile();
        //this.creatUltraResPAVFPGMap();
        //this.creatPAVPosFile();
        //this.creatPAVUltraresPosFile();
        
    }
    
    public void jointLinkagePipe () {
//*************PAV identification section******************
        //this.mergeTBTByRow();
        
        //this.sliceTBT();
        
        //this.creatTagPMap();
        
        //this.mapTBTJoint(); //impossible in this computer
        
        //this.mergeMappingResult(true);
        
        //this.creatTagJointPGMap();
        
        //this.creatPredictFile(3);
        
        this.updateResolution(3);
        
        
        
        //this.checkJointMappingQuality();
        
        
        //this.filteredTagJointPGMapByThresh(); //Deprecated
        
        //this.jointStatistics(); //Deprecated
        
        //this.alignTagOfJointPGMap(); //Deprecated
        
        //this.identifyPAV(); //Deprecated
        
        
//*******************Standerdize data**********************
        //this.getB73TagPosCount();
    }
    
    public void assoPipe () {//LRatio > 4, 90% on chromosome, 60% 100kb region, population structure affect the mapping
        //this.mergeTBTByRow();
        
        //this.sliceTBT();
        
        //this.mapTBTGWAS(); //impossible in this computer
        
        //this.mergeMappingResult(false);
        
        //this.creatTagGWASPGMap();
        
        //this.checkGWASMappingQuality();
        
        //this.mergePGMaps();
        
        //this.creatPredictFile(1);
        
        //this.updateResolution(1);
    }
    
    public void figurePipe () {
        //this.figureGMappingQuality();
        //this.figureJMappingQuality();
        //this.figureFPGMapQuality();
        //this.figureResolution();
        //this.fitureFPGResolution();
        this.enrichmentLowCopyNum();
        //this.enrichmentGeneVsRepeat();
        //this.readCountAlongChr();
        //this.deletionCNVMedianRatio(); //Deprecated
    }
    
    public void trainPipe () {
        this.creatUB73TBT();
        //this.sliceUB73TBT();
        //this.mapUB73TBTGWAS();
        //this.mergeUB73MappingResult();
        //this.creatUB73TagGWASPGMap();
        //this.creatUB73TagMPGMap();
        //this.makeTrainingData();
    }
    
    public void callAmes282GenotypeVariableBin () {
        double pavP = 0.001; double cnvP = 0.001; double occurence = 0.01;
        String refName = "b73";
        String VBBTFileS = new File (fileOfData[15], "/vbbt/ames282.ataxa.vbbt.50.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[15], "/hapmap/ames282.ataxa.50.p001.01.pav.hapmap.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[15], "/hapmap/ames282.ataxa.50.p001.01.cnv.hapmap.txt").getAbsolutePath();
        VariableBinByTaxa vbbt = new VariableBinByTaxa (VBBTFileS);
        RatioMatrixVariableBin rmvb = new RatioMatrixVariableBin (vbbt, refName, pavP, cnvP, occurence);
        rmvb.mkPAVgenotypeByFDR(pavHapMap, 0.3);
        rmvb.mkCNVgenotypeByFDR(cnvHapMap, 0.3);
        //rmvb.mkMedianRatioFile(cnvHapMap);
        System.out.println(rmvb.getPAVPortion());
        System.out.println(rmvb.getCNVPortion());
        System.out.println(rmvb.getPAVOrCNVPortion());
    }
    
    public void call282GenotypeVariableBin () {
        double pavP = 0.001; double cnvP = 0.001; double occurence = 0.01;
        String refName = "b73";
        String VBBTFileS = new File (fileOfData[15], "/vbbt/282.merge.vbbt.50.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[15], "/hapmap/282.merge.50.p001.05.pav.hapmap.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[15], "/hapmap/282.merge.50.p001.05.cnv.hapmap.txt").getAbsolutePath();
        VariableBinByTaxa vbbt = new VariableBinByTaxa (VBBTFileS);
        RatioMatrixVariableBin rmvb = new RatioMatrixVariableBin (vbbt, refName, pavP, cnvP, occurence);
        //rmvb.mkPAVgenotypeByFDR(pavHapMap, 0.15);
        rmvb.mkPAVgenotypeByP(pavHapMap, pavP);
        rmvb.mkCNVgenotypeByFDR(cnvHapMap, 0.3);
        String medianRatioFileS = new File (fileOfData[15], "/hapmap/paper.medianRatio.txt").getAbsolutePath();
        //rmvb.mkMedianRatioFile(medianRatioFileS);
        
        System.out.println(rmvb.getPAVPortion());
        System.out.println(rmvb.getCNVPortion());
        System.out.println(rmvb.getPAVOrCNVPortion());
    }
    
    public void callPaperGenotypeVariableBin () {
        double pavP = 0.01; double cnvP = 0.01; double occurence = 0.03;
        String refName = "b73";
        String VBBTFileS = new File (fileOfData[15], "/vbbt/paper.merge.vbbt.50.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[15], "/hapmap/paper.merge.50.p001.03.pav.hapmap.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[15], "/hapmap/paper.merge.50.p001.03.cnv.hapmap.txt").getAbsolutePath();
        VariableBinByTaxa vbbt = new VariableBinByTaxa (VBBTFileS);
        RatioMatrixVariableBin rmvb = new RatioMatrixVariableBin (vbbt, refName, pavP, cnvP, occurence);
        //rmvb.mkPAVgenotypeByFDR(pavHapMap, 0.05);
        //rmvb.mkCNVgenotypeByFDR(cnvHapMap, 0.05);
        rmvb.mkPAVgenotypeByP(pavHapMap, pavP);
        rmvb.mkCNVgenotypeByP(cnvHapMap, cnvP);
        String medianRatioFileS = new File (fileOfData[15], "/hapmap/paper.medianRatio.txt").getAbsolutePath();
        //rmvb.mkMedianRatioFile(medianRatioFileS);
        
        System.out.println(rmvb.getPAVPortion());
        System.out.println(rmvb.getCNVPortion());
        System.out.println(rmvb.getPAVOrCNVPortion());
    }
    
    public void callNAMGenotypeVariableBin () {
        double pavP = 0.001; double cnvP = 0.001; double occurence = 0.2;
        String refName = "b73";
        String VBBTFileS = new File (fileOfData[15], "/vbbt/NamParents.merge.vbbt.60.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[15], "/hapmap/namParents.merge.60.p001.01.pav.hapmap.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[15], "/hapmap/namParents.merge.60.p001.01.cnv.hapmap.txt").getAbsolutePath();
        VariableBinByTaxa vbbt = new VariableBinByTaxa (VBBTFileS);
        RatioMatrixVariableBin rmvb = new RatioMatrixVariableBin (vbbt, refName, pavP, cnvP, occurence);
        rmvb.mkPAVgenotypeByFDR(pavHapMap, 0.3);
        rmvb.mkCNVgenotypeByFDR(cnvHapMap, 0.3);
        
        String ratioFileS = new File (fileOfData[15], "/hapmap/namParents.ratio.txt").getAbsolutePath();
        rmvb.mkRatioFile(ratioFileS);
        System.out.println(rmvb.getPAVPortion());
        System.out.println(rmvb.getCNVPortion());
        System.out.println(rmvb.getPAVOrCNVPortion());
        
        //**************************
        //String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
        //String smoothRatioFileS = new File (fileOfData[15], "/hapmap/namParents.smoothratio.txt").getAbsolutePath();
        //SmoothRatioMatrix srm = new SmoothRatioMatrix(genomeInfoFileS, 1000000, rmvb);
        //srm.writeRatioFileS(smoothRatioFileS);
        //**************************
    }
    
    public void selectVariableBBT () {
        VariableBinByTaxa bbt;
        
        String paperTaxaFileS = new File (fileOfData[15], "/taxaName/paper.all.merge.lower.taxa.txt").getAbsolutePath();
        String binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/all.merge.vbbt.50.bin").getAbsolutePath();
        String paperVBBTFileS = new File (fileOfData[15], "/vbbt/paper.merge.vbbt.50.bin").getAbsolutePath();
        //bbt = new VariableBinByTaxa (binByTaxaAb73FileS);
        //bbt.selectTaxa(paperTaxaFileS);
        //bbt.writeTxtFile(paperVBBTFileS);
        //bbt.writeBinaryFile(paperVBBTFileS);
        
        String ames282TaxaFileS = new File (fileOfData[15], "/taxaName/282Ames.txt").getAbsolutePath();
        binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/all.ataxa.vbbt.50.bin").getAbsolutePath();
        String ames282VBBTFileS = new File (fileOfData[15], "/vbbt/ame282.ataxa.vbbt.50.bin").getAbsolutePath();
        //bbt = new VariableBinByTaxa (binByTaxaAb73FileS);
        //bbt.selectTaxa(ames282TaxaFileS);
        //bbt.writeTxtFile(ames282VBBTFileS);
        //bbt.writeBinaryFile(ames282VBBTFileS);
        
        String namParentsFileS = new File (fileOfData[15], "/taxaName/NAMparents.txt").getAbsolutePath();
        String binByTaxaMergeFileS = new File (fileOfData[15], "/vbbt/all.merge.vbbt.50.bin").getAbsolutePath();
        String namParentsVBBTFileS = new File (fileOfData[15], "/vbbt/namParents.merge.vbbt.50.bin").getAbsolutePath();
        //bbt = new VariableBinByTaxa (binByTaxaMergeFileS);
        //bbt.selectTaxa(namParentsFileS);
        //bbt.writeTxtFile(namParentsVBBTFileS);
        //bbt.writeBinaryFile(namParentsVBBTFileS);
        
        String assoTaxaFileS = new File (fileOfData[15], "/taxaName/282.txt").getAbsolutePath();
        binByTaxaMergeFileS = new File (fileOfData[15], "/vbbt/all.merge.vbbt.50.bin").getAbsolutePath();
        String assoVBBTFileS = new File (fileOfData[15], "/vbbt/282.merge.vbbt.50.bin").getAbsolutePath();
        bbt = new VariableBinByTaxa (binByTaxaMergeFileS);
        bbt.selectTaxa(assoTaxaFileS);
        //bbt.writeTxtFile(assoVBBTFileS);
        bbt.writeBinaryFile(assoVBBTFileS);
    }
    
    public void mergeTaxaVariableBBTNAM () {//merge multiple samples
        String binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/NamParents.ataxa.vbbt.60.bin").getAbsolutePath();
        String binByTaxaMergeFileS = new File (fileOfData[15], "/vbbt/NamParents.merge.vbbt.60.bin").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (binByTaxaAb73FileS);
        bbt.mergeByTaxaName();
        bbt.writeBinaryFile(binByTaxaMergeFileS);
        bbt.screePrintAverCount();
    }
    
    public void mergeTaxaVariableBBT () {//merge multiple samples
        String binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/all.ataxa.vbbt.50.bin").getAbsolutePath();
        String binByTaxaMergeFileS = new File (fileOfData[15], "/vbbt/all.merge.vbbt.50.bin").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (binByTaxaAb73FileS);
        bbt.mergeByTaxaName();
        bbt.writeBinaryFile(binByTaxaMergeFileS);
        bbt.screePrintAverCount();
    }
    
    public void addRemoveb73VairableBBTNAM () {//add "b73" taxa, remove other B73s
        String binByTaxaFileS = new File (fileOfData[15], "/vbbt/NamParents.vbbt.60.bin").getAbsolutePath();
        String binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/NamParents.ataxa.vbbt.60.bin").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (binByTaxaFileS);
        bbt.addb73Taxon();
        bbt.writeBinaryFile(binByTaxaAb73FileS);
        bbt.screePrintAverCount();
    }
    
    public void addRemoveb73VairableBBT () {//add "b73" taxa, remove other B73s
        String binByTaxaFileS = new File (fileOfData[15], "/vbbt/all.vbbt.50.bin").getAbsolutePath();
        String binByTaxaAb73FileS = new File (fileOfData[15], "/vbbt/all.ataxa.vbbt.50.bin").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (binByTaxaFileS);
        bbt.addb73Taxon();
        bbt.writeBinaryFile(binByTaxaAb73FileS);
        bbt.screePrintAverCount();
    }
    
    public void taxaStatistics () {
        
        String paperTaxaNameFileS = new File (fileOfData[15], "/taxaName/paper.full.taxa.txt").getAbsolutePath();
        String mergeTaxaNameFileS = new File (fileOfData[15], "/taxaName/paper.merge.taxa.txt").getAbsolutePath();
        String lowerMergeTaxaNameFileS = new File (fileOfData[15], "/taxaName/paper.merge.lower.taxa.txt").getAbsolutePath();
        String lower2FullTaxaNameFileS = new File (fileOfData[15], "/taxaName/paper.merge.lower2Full.taxa.txt").getAbsolutePath();
        TaxaNameUtils2 tnu = new TaxaNameUtils2 (paperTaxaNameFileS);
        tnu.mkMergeTaxaNameFileS(mergeTaxaNameFileS, lowerMergeTaxaNameFileS); //multiple samples in one taxa
        tnu.mkLower2FullNameFileS(lowerMergeTaxaNameFileS, lower2FullTaxaNameFileS);
        
        String paperKeyFileS = new File (fileOfData[15], "/taxaName/source/paper_key.txt").getAbsolutePath();
        tnu = new TaxaNameUtils2 (mergeTaxaNameFileS);
        String panelTaxaNumTableS = new File (fileOfData[15], "/taxaName/panelTaxaNum.txt").getAbsolutePath();
        tnu.mkPanelTaxaNumTable(paperKeyFileS, panelTaxaNumTableS);
    }
    
    public void selectTaxaPanel () {
        //working on NAM, Ames, BREAD, CN-NAM, 282, IBM, Imputation p1
        String binByTaxaFileS = new File (fileOfData[15], "/vbbt/all.vbbt.50.bin").getAbsolutePath();
        String allTaxaName = new File (fileOfData[15], "/taxaName/Jan.taxa.fullname.txt").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (binByTaxaFileS);
        bbt.mkTaxaNameTable(allTaxaName);
        String paperKeyFileS = new File (fileOfData[15], "/taxaName/source/paper_key.txt").getAbsolutePath();
        String paperTaxaNameFileS = new File (fileOfData[15], "/taxaName/paper.full.taxa.txt").getAbsolutePath();
        bbt.mkPaperTaxaNameTable(paperKeyFileS, paperTaxaNameFileS); //find shared taxa between keyfile and Janurary build
    }
    public void creatVariableBinByTaxaNAM () {
        String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String mergedTBTFileS = "M:/pav/mergedTBT/tbt_NamParents.byte";
        String binSizeFileS = new File (fileOfData[15], "/vbbt/binSize_NamParents.txt").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (genomeInfoFileS, tagPMapFileS, binSizeFileS, mergedTBTFileS);
        String binByTaxaFileS = new File (fileOfData[15], "/vbbt/NamParents.vbbt.60.bin").getAbsolutePath();
        bbt.writeBinaryFile(binByTaxaFileS);
    }
    
    public void creatVariableBinByTaxa () {
        String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String mergedTBTFileS = "J:/tbt_20.byte";
        //String mergedTBTFileS = "M:/pav/mergedTBT/tbt_test.byte";
        String binSizeFileS = new File (fileOfData[15], "/vbbt/binSize_50.txt").getAbsolutePath();
        VariableBinByTaxa bbt = new VariableBinByTaxa (genomeInfoFileS, tagPMapFileS, binSizeFileS, mergedTBTFileS);
        
        String binByTaxaFileS = new File (fileOfData[15], "/vbbt/all.vbbt.50.bin").getAbsolutePath();
        bbt.writeBinaryFile(binByTaxaFileS);
    }
    
    public void creatBinSize () {
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String binSizeFileS = new File (fileOfData[15], "/vbbt/binSize.txt").getAbsolutePath();
        int taxaNum = 22743;
        double averageCount = 50;
        VariableBinByTaxa vbbt = new VariableBinByTaxa ();
        vbbt.mkBinSizeFile(tagCountFileS, tagPMapFileS, binSizeFileS, taxaNum, averageCount);
    }
    
    public void creatBinSizeNAMParents () {
        String allTaxaName = new File (fileOfData[15], "/taxaName/all.taxa.txt").getAbsolutePath();
        String namParentFileS = new File (fileOfData[15], "/taxaName/NAMparents.txt").getAbsolutePath();
        String namParentFullNameFileS = new File (fileOfData[15], "/taxaName/NAMParents_full.txt").getAbsolutePath();
        TaxaNameUtils tnu = new TaxaNameUtils (allTaxaName);
        //tnu.mkNAMParentsFullName(namParentFileS, namParentFullNameFileS);
        
        //String inputTBTFileS = "J:/tbt_20.byte";
        
        String outputTBTFileS = "M:/pav/mergedTBT/tbt_NamParents.byte";
        PAVUtils util = new PAVUtils ();
        //util.mkSubTBTByTaxa(inputTBTFileS, namParentFullNameFileS, outputTBTFileS);
        
        String namParentTagCountFileS = new File (fileOfData[2], "/NamParents_full.cnt").getAbsolutePath();
        //util.convertTBTToTagCount(outputTBTFileS, namParentTagCountFileS);
        
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String binSizeFileS = new File (fileOfData[15], "/vbbt/binSize_NamParents.txt").getAbsolutePath();
        VariableBinByTaxa vbbt = new VariableBinByTaxa ();
        vbbt.mkBinSizeFile(namParentTagCountFileS, tagPMapFileS, binSizeFileS, 27, 60);
    }
    
    public void figureFPGMapQuality () {
        String finalPGFileS = new File (fileOfData[8], "All.fpgmap.res.txt").getAbsolutePath();
        String outfileS = new File (fileOfData[11], "quality/quality_fpg.txt").getAbsolutePath();
        TagFPGMap fpgmap = new TagFPGMap (finalPGFileS);
        int outNum = 20000;
        fpgmap.mkFPGMappingQuality(outfileS, outNum);
    }
    
    public void figureJMappingQuality () {
        double pThresh = 0.01;
        int outNum = 20000;
        String name = "quality/quality_J_p"+String.valueOf(pThresh)+".txt";
        String tagJointPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String outfileS = new File (fileOfData[11], name).getAbsolutePath();
        TagJointPGMap jpgmap = new TagJointPGMap(tagJointPGMapFileS, false);
        jpgmap.mkJMappingQuality(outfileS, pThresh, outNum);
    }
    
    public void figureGMappingQuality () {
        double pTresh = 1E-2, lThresh = 0;
        String name = "quality/quality_G_l"+String.valueOf(lThresh)+".txt";
        String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
        String outfileS = new File (fileOfData[11], name).getAbsolutePath();
        
        int outNum = 20000;
        TagGWASPGMap gpgmap = new TagGWASPGMap(tagGWASPGMapFileS);
        gpgmap.mkGMappingQuality(outfileS, pTresh, lThresh, outNum);
    }
    
    public void deletionCNVMedianRatio () {
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/NAMparents.50k.p001.o03.rmx.bin").getAbsolutePath();
        String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/282Ames.200k.p001.o03.rmx.bin").getAbsolutePath();
        String medianRatioFileS = "E:/Research/pav/attDis/delCnvMedian.txt";
        RatioMatrix rm = new RatioMatrix (pMatrixFileS);
        rm.mkMedianRatioFile(medianRatioFileS);
    }
    
    public void readCountAlongChr () {
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String readCountChrFileS = "E:/Research/pav/binSize/readCountChr.txt";
        int chr = 1;
        PAVUtils util = new PAVUtils();
        util.mkReadCountOnChrFile(tagCountFileS, tagPMapFileS, readCountChrFileS, chr);
    }
    
    public void enrichmentGeneVsRepeat () {
        PAVUtils util = new PAVUtils();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String hitNum26MFileS = new File (fileOfData[11], "enrichment/26MB73_hitNum.txt").getAbsolutePath();
        String fastaP26MRefFileS = new File (fileOfData[11], "enrichment/p26MB73.fasta").getAbsolutePath();
        util.mkRefHitNumFastaFile(tagCountFileS, hitNum26MFileS, fastaP26MRefFileS);
        String hitNumFPGFileS = new File (fileOfData[11], "enrichment/fpgB73_hitNum.txt").getAbsolutePath();
        String fastaPFPGRefFileS = new File (fileOfData[11], "enrichment/pFPGB73.fasta").getAbsolutePath();
        util.mkRefHitNumFastaFile(tagCountFileS, hitNumFPGFileS, fastaPFPGRefFileS);
        
    }
    
    public void enrichmentLowCopyNum () {
        String finalPGMapFileS = new File (fileOfData[8], "All.fpgmap.res.txt").getAbsolutePath();
        //TagFPGMap fpg = new TagFPGMap(finalPGMapFileS);
        String samFileS = new File (fileOfData[6], "fpg-very-sensitive-local-k20.sam").getAbsolutePath();
        //TagPMap pmap = new TagPMap(fpg, samFileS);
        String tagFPGPMapFileS = new File (fileOfData[6], "fpg-k20.pmap").getAbsolutePath();
        //pmap.writeTagPMap(tagFPGPMapFileS);
        
        
        int tagLimit = 200000;
        PAVUtils util = new PAVUtils();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        
        String tag26MPMapFileS = new File (fileOfData[6], "26Mtag-k20.pmap").getAbsolutePath();
        String hitNum26MFileS = new File (fileOfData[11], "enrichment/26MB73_hitNum.txt").getAbsolutePath();
        //util.mkRefPerfectHitNumFile(tag26MPMapFileS, tagCountFileS, hitNum26MFileS, tagLimit);
        
        
        String hitNumFPGFileS = new File (fileOfData[11], "enrichment/fpgB73_hitNum.txt").getAbsolutePath();
        util.mkRefPerfectHitNumFile(tagFPGPMapFileS, tagCountFileS, hitNumFPGFileS, tagLimit);
    }
    
    public void fitureFPGResolution () {
        String finalPGFileS = new File (fileOfData[8], "All.fpgmap.res.txt").getAbsolutePath();
        String resolutionFileS = new File (fileOfData[11], "resolution/resolution_FPG.txt").getAbsolutePath();
        TagFPGMap fpgmap = new TagFPGMap (finalPGFileS);
        int outnum = 20000;
        fpgmap.mkResolutionFile(outnum, resolutionFileS, false);
    }
    
    public void figureResolution () {
        String UB73MPGMapFileS = new File (fileOfData[8], "UB73.mpgmap.txt").getAbsolutePath();
        String resolutionFileS = new File (fileOfData[11], "resolution/resolution_GJ.txt").getAbsolutePath();
        TagMPGMap mpgmap = new TagMPGMap(UB73MPGMapFileS);
        int outnum = 20000;
        double lCutoff = 4, pCutoff = 0.01;
        mpgmap.mkResolutionFile(lCutoff, pCutoff, outnum, resolutionFileS);
    }
    
    public void makeTrainingData () {
        String tagMPGMapFileS = new File (fileOfData[8], "UB73.mpgmap.txt").getAbsolutePath();
        String B73TrainingPGMapFileS = new File (fileOfData[10], "UB73train.mpgmap.txt").getAbsolutePath();
        TagMPGMap mpgmap = new TagMPGMap (tagMPGMapFileS);
        mpgmap.writeTxtFile(B73TrainingPGMapFileS, 2000, 30000);
        
        
        //output training data
        String tagHitFileS = new File (fileOfData[10], "tag_hits_k10_very_sensi.txt").getAbsolutePath();
        String recombinationFileS = new File (fileOfData[10], "recombination.txt").getAbsolutePath();
        String trainingFileS = new File (fileOfData[10], "UB73train.csv").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        Train tr = new Train (B73TrainingPGMapFileS);
        tr.addTagCountHits(tagCountFileS, tagHitFileS);
        tr.addRecombination(recombinationFileS);
        tr.preTransform();
        tr.logTransform();
        //tr.rankTransform();
        //tr.boxcoxTransform();
        tr.writeTrainingFile(trainingFileS);
    }
    
    public void creatUB73TagMPGMap () {
        String tagGWASPGMapFileS = new File (fileOfData[8], "UB73.gpgmap.txt").getAbsolutePath();
        String tagJointPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String tagMPGMapFileS = new File (fileOfData[8], "UB73.mpgmap.txt").getAbsolutePath();
        TagMPGMap mpgmap = new TagMPGMap (tagGWASPGMapFileS, tagJointPGMapFileS);
        mpgmap.writeTxtFile(tagMPGMapFileS);
    }
    
    public void creatUB73TagGWASPGMap () {
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String tagGWASGMapFileS = fileOfData[5].getAbsolutePath() + "/mergedUB73.mapping.gwas.txt";
        String tagGWASPGMapFileS = new File (fileOfData[8], "UB73.gpgmap.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        TagGWASPGMap gpgmap = new TagGWASPGMap (tagPMapFileS, tagGWASGMapFileS, tagCountFileS);
        gpgmap.writeTxtFile(tagGWASPGMapFileS);
    }
    
    public void mergeUB73MappingResult () {
        String inputDirS = new File (fileOfData[5], "UB73").getAbsolutePath();
        String mergedFileS = fileOfData[5].getAbsolutePath() + "/mergedUB73.mapping.gwas.txt";
        new PAVUtils().mergeGeneticMappingFiles(inputDirS, mergedFileS);
    }
    
    public void mapUB73TBTGWAS () {
        String tbtDirS = new File (fileOfData[1], "sample").getAbsolutePath();
        String anchorMapFileS = new File (fileOfData[4], "testData/impAll/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String mappingResultDirS = "M:/test";
        
        TagAgainstAnchorLongTimePosBlock taa = new TagAgainstAnchorLongTimePosBlock (tbtDirS, anchorMapFileS, tagPMapFileS, mappingResultDirS);
    }
    
    public void sliceUB73TBT () {
        String UB73TBT = new File (fileOfData[1], "UB73.tbt.byte").getAbsolutePath();
        String sliceTBTDirS = new File (fileOfData[3], "UB73").getAbsolutePath();
        int sliceNum = 300;
        PAVUtils uti = new PAVUtils();
        uti.sliceTBT(sliceTBTDirS, UB73TBT, sliceNum);
    }
    
    public void creatUB73TBT () {
        String mergedTBT = "I:/tbt_20.byte";
        String UB73TBT = new File (fileOfData[1], "UB73.tbt.byte").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        PAVUtils util = new PAVUtils();
        util.mkURefTagTBTFile(tagPMapFileS, tagCountFileS, mergedTBT, UB73TBT);
    }
    
    public void calCnvRatioByB73 () {
        String cnvPosFileS = new File (fileOfData[9], "/tagRatio/cnv_posi.txt").getAbsolutePath();
        String cnvRatioFileS = new File (fileOfData[9], "/tagRatio/cnv_ratio.txt").getAbsolutePath();
        
    }
    
    public void assignPosToCnvTBTAndMergeTaxa () {//just a test, this method consumes too much memory, can't sort
        String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
        String cnvPosFileS = new File (fileOfData[9], "/tagPos/cnv_pos.txt").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String totalCountFileS = new File (fileOfData[9], "/totalCount/taxaCount.txt").getAbsolutePath();
        PAVUtils util = new PAVUtils();
        util.makeCNVPosFile(cnvTBTFileS, cnvPosFileS, tagPMapFileS, totalCountFileS);
        
    }
    
    public void getTaxaTagCountInBin () {
        String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
        String binaryBinFileS = new File (fileOfData[9], "/tagCountInBin/binCount.bin").getAbsolutePath();
        TaxaTagCountInBin tagBin = new TaxaTagCountInBin(genomeInfoFileS,cnvTBTFileS, tagPMapFileS, 50000);
        tagBin.writeBinaryFile(binaryBinFileS);
        
        String txtBinFileS = new File (fileOfData[9], "/tagCountInBin/binCount.txt").getAbsolutePath();
        TaxaTagCountInBin tagBin2 = new TaxaTagCountInBin(binaryBinFileS);
        tagBin2.mergeTaxa();
        tagBin2.writeTxtFile(txtBinFileS);
    }
    
    public void getTBTOfCNVTag () {
        String cnvTagCountFileS =  new File (fileOfData[9], "/tagCount/cnv.cnt").getAbsolutePath();
        String mergedTBTFileS = "Z:/tbt_20.byte";
        //String mergedTBTFileS = new File (fileOfData[1], "tbt_test.byte").getAbsolutePath();
        String taxaNameFileS = new File (fileOfData[9], "/taxaName/interet_sample.txt").getAbsolutePath();
        String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv.tbt.byte").getAbsolutePath();
        PAVUtils util = new PAVUtils();
        util.makeSubTBTByTaxaTagCount(mergedTBTFileS, cnvTagCountFileS, taxaNameFileS, cnvTBTFileS);
    }
    
    
    public void mkPAVgenotypeFrequency () {
        String pavHapMap = new File (fileOfData[9], "/genotype/NAM12.hmp.txt").getAbsolutePath();
        String frequencyFileS = "E:/Research/pav/pavGenotype/cnv_pav.NAM12.fre.txt";
        HapMapUtils util = new HapMapUtils ();
        util.mkAlleleFrequency(pavHapMap, frequencyFileS);
    }
    
    public void callCNVgenotypes () {
        
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/282Ames.200k.p001.o03.rmx.bin").getAbsolutePath();
        //String cnvHapMap = new File (fileOfData[9], "/genotype/282Ames.cnv.200k.p001.o03.hmp.txt").getAbsolutePath();
        
        String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/paper.high.50k.unknown.rmx.bin").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[9], "/genotype/paper.high.cnv.50k.unknown.hmp.txt").getAbsolutePath();
        
        
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/NAMparents.50k.p001.o03.rmx.bin").getAbsolutePath();
        //String cnvHapMap = new File (fileOfData[9], "/genotype/NAMparents.cnv.50k.p001.o03.hmp.txt").getAbsolutePath();
        
        RatioMatrix pm = new RatioMatrix(pMatrixFileS);
        System.out.println(pm.getPAVLociAll());
        System.out.println(pm.getCNVLociAll());
        System.out.println(pm.getPAVOrCNVLociAll());
        
        pm.mkCNVgenotypeByP(cnvHapMap, 0.3);
    }
    
    public void callPAVGenotypes () {
        
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/NAMparents.50k.p01.o03.rmx.bin").getAbsolutePath();
        //String pavHapMap = new File (fileOfData[9], "/genotype/NAMparents.pav.50k.p01.o03.hmp.txt").getAbsolutePath();
        
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/282Ames.200k.p005.o03.rmx.bin").getAbsolutePath();
        //String pavHapMap = new File (fileOfData[9], "/genotype/282Ames.pav.200k.p005.o03.hmp.txt").getAbsolutePath();
        
        String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/paper.high.50k.unknown.rmx.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[9], "/genotype/paper.high.pav.50k.unknown.hmp.txt").getAbsolutePath();
        
        //String pMatrixFileS = new File (fileOfData[9], "/ratioMatrix/282Ames.200k.pcut.rmx.bin").getAbsolutePath();
        //String pavHapMap = new File (fileOfData[9], "/genotype/282Ames.pav.200k.pcut.hmp.txt").getAbsolutePath();
        
        //String pavMedianRatioFileS = new File (fileOfData[9], "/genotype/282Ames.500k.p001.o01.medianratio.txt").getAbsolutePath();
        //String pavAverRatioFileS = new File (fileOfData[9], "/genotype/282Ames.200k.p001.o03.averageratio.txt").getAbsolutePath();
        
        RatioMatrix pm = new RatioMatrix(pMatrixFileS);
        pm.mkPAVgenotypeByP(pavHapMap, 0.3);
        
        String svCountFileS = new File (fileOfData[9], "/ratioMatrix/svCount.txt").getAbsolutePath();
        pm.mkSVCountInRegion(1000000, svCountFileS);
        
        System.out.println(pm.getPAVLociAll());
        System.out.println(pm.getCNVLociAll());
        System.out.println(pm.getPAVOrCNVLociAll());
        
        //pm.mkPAVgenotypeMedianRatioByP(pavMedianRatioFileS, 0.3);
        //pm.mkPAVgenotypeAverageRatioByP(pavAverRatioFileS, 0.3);
    }
    
    public void creatRatioMatrix () {
        
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.NAMparents.mtaxa.50k.bbt.txt").getAbsolutePath();
        String ratioFileS = new File (fileOfData[9], "/ratioMatrix/NAMparents.50k.p001.o03.rmx.bin").getAbsolutePath();
        
        //String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.282Ames.ataxa.200k.bbt.txt").getAbsolutePath();
        //String ratioFileS = new File (fileOfData[9], "/ratioMatrix/282Ames.200k.pcut.rmx.bin").getAbsolutePath();
        
        //String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/paper.high.mtaxa.50k.bbt.txt").getAbsolutePath();
        //String ratioFileS = new File (fileOfData[9], "/ratioMatrix/paper.high.50k.unknown.rmx.bin").getAbsolutePath();
        
        
        BinByTaxa bbt = new BinByTaxa(binByTaxaMTaxaFileS);
        //ouble pavP = 0.0048, cnvP = 0.000137, occurence = 0.01; //fdr = 0.01. 200k;
        //double pavP = 0.0029, cnvP = 0.000424, occurence = 0.03; //fdr = 0.01. 500k;
        double pavP = 0.001, cnvP = 0.001, occurence = 0.03;
        //RatioMatrix pm = new RatioMatrix(bbt, "b73", fdr, occurence);
        RatioMatrix pm = new RatioMatrix(bbt, "b73", pavP, cnvP, occurence);
        //pm.displayFDRpCutoff(0.1);
        pm.writeBinaryFile(ratioFileS);
        
        System.out.println(pm.getPAVLociAll());
        System.out.println(pm.getCNVLociAll());
        System.out.println(pm.getPAVOrCNVLociAll());
    }
    
    
    public void creatRatioMatrixByZscore () {
        String zscoreTableS = new File (fileOfData[9], "/ratioMatrix2/cutoff/zscore.200k.txt").getAbsolutePath();
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/paper.mtaxa.200k.bbt.txt").getAbsolutePath();
        String ratioFileS = new File (fileOfData[9], "/ratioMatrix/paper.200k.zscore.rmx.bin").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaMTaxaFileS);
        double zCutoff = 5, occurence = 0.03;
        RatioMatrix pm = new RatioMatrix(bbt, "b73", zscoreTableS, zCutoff, occurence);
        System.out.println(pm.getPAVLociAll());
        System.out.println(pm.getCNVLociAll());
        System.out.println(pm.getPAVOrCNVLociAll());
        pm.writeBinaryFile(ratioFileS);
    }
    
    public void creatGenotypeByRDis () {
        String refName = "b73";
        String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.282Ames.ataxa.500k.bbt.txt").getAbsolutePath();
        String pavPCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/pav.cutoff.500k.txt").getAbsolutePath();
        String cnvPCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/cnv.cutoff.500k.txt").getAbsolutePath();
        String pavHapMap = new File (fileOfData[9], "/genotype2/282Ames.pav.500k.pcut.hmp.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[9], "/genotype2/282Ames.cnv.500k.pcut.hmp.txt").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaSelFileS);
        double occurence = 0.05;
        RatioMatrix rm = new RatioMatrix (bbt, refName, pavPCutoffFileS, cnvPCutoffFileS, occurence, false);
        
        rm.mkPAVgenotypeByP(pavHapMap, 0.3);
        rm.mkCNVgenotypeByP(cnvHapMap, 0.3);
        
        System.out.println(rm.getPAVLociAll());
        System.out.println(rm.getCNVLociAll());
        System.out.println(rm.getPAVOrCNVLociAll());
        
    }
    
    public void creatGenotypeByPDis () {
        String ratioFileS = new File (fileOfData[9], "/ratioMatrix2/282Ames.500k.pcut.rmx.bin").getAbsolutePath();
        String pavHapMap = new File (fileOfData[9], "/genotype2/282Ames.pav.500k.pcut.hmp.txt").getAbsolutePath();
        String cnvHapMap = new File (fileOfData[9], "/genotype2/282Ames.cnv.500k.pcut.hmp.txt").getAbsolutePath();
        RatioMatrix pm = new RatioMatrix(ratioFileS);
        pm.mkPAVgenotypeByP(pavHapMap, 0.3);
        pm.mkCNVgenotypeByP(cnvHapMap, 0.3);
        
        System.out.println(pm.getPAVLociAll());
        System.out.println(pm.getCNVLociAll());
        System.out.println(pm.getPAVOrCNVLociAll());
    }
    
    public void creatRatioMatrixByPDis () {
        String refName = "b73";
        String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.282Ames.ataxa.500k.bbt.txt").getAbsolutePath();
        String pavPCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/pav.cutoff.500k.txt").getAbsolutePath();
        String cnvPCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/cnv.cutoff.500k.txt").getAbsolutePath();
        String ratioFileS = new File (fileOfData[9], "/ratioMatrix2/282Ames.500k.pcut.rmx.bin").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaSelFileS);
        double occurence = 0.05;
        RatioMatrix rm = new RatioMatrix (bbt, refName, pavPCutoffFileS, cnvPCutoffFileS, occurence, true);
        rm.writeBinaryFile(ratioFileS);
        
        System.out.println(rm.getPAVLociAll());
        System.out.println(rm.getCNVLociAll());
        System.out.println(rm.getPAVOrCNVLociAll());
    }
    
    public void permutationCutoff () {
        String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.B73.ataxa.200k.bbt.txt").getAbsolutePath(); //one taxa is all 0, needs to be removed
        BinByTaxa bbt = new BinByTaxa(binByTaxaSelFileS);
        RatioMatrixB73Permute pm = new RatioMatrixB73Permute (bbt, "b73");
        //double fdr = 0.01;
        //double[] pCut = pm.getPAtFDR(fdr);
        //System.out.println("PAV ratio cutoff is " + pCut[0]);
        //System.out.println("CNV ratio cutoff is " + pCut[1]);
        
        
        //*********************************************************************
        String pavCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/pav.cutoff.200k.txt").getAbsolutePath();
        String cnvCutoffFileS = new File (fileOfData[9], "/ratioMatrix2/cutoff/cnv.cutoff.200k.txt").getAbsolutePath();
        
        //pm.mkPAVCutoffTable(pavCutoffFileS, fdr);
        //pm.mkCNVCutoffTable(cnvCutoffFileS, fdr);
        
        int sampleNum = 1, popSize = 1000;
        String zscoreTableS = new File (fileOfData[9], "/ratioMatrix2/cutoff/zscore.200k.txt").getAbsolutePath();
        
        //RatioMatrixB73Permute b73pm = pm.getPermutedB73Ratiomatrix(bbt, sampleNum, popSize);
        pm.mkZscoreTable(zscoreTableS);
        //*********************************************************************
    }
    
    public void selectTaxaInBBT () {
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        String taxaNameFileS = new File (fileOfData[9], "/taxaName/paper.high.merge.lower.taxa.txt").getAbsolutePath();
        String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/paper.high.mtaxa.200k.txt").getAbsolutePath();
        
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/NAM.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.NAM.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/282.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.282.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/Ames.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.Ames.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/NAMparents.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.NAMparents.mtaxa.bbt.txt").getAbsolutePath();
        
        //String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.ataxa.200k.bbt.txt").getAbsolutePath();
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/282Ames.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.282Ames.ataxa.200k.bbt.txt").getAbsolutePath();
        
        //String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.ataxa.200k.bbt.txt").getAbsolutePath();
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/B73.txt").getAbsolutePath();
        //String binByTaxaSelFileS = new File (fileOfData[9], "/bbt/cnv.B73.ataxa.200k.bbt.txt").getAbsolutePath();
        
        BinByTaxa bbt = new BinByTaxa(binByTaxaMTaxaFileS);
        System.out.println(bbt.getTagCountAll());
        bbt.selectTaxa(taxaNameFileS);
        bbt.writeTxtFile(binByTaxaSelFileS);
    }
    
    public void creatTaxaNameList () {
        //String taxaNameFileS = new File (fileOfData[9], "/taxaName/mergeTaxa.txt").getAbsolutePath();
        String panel282FileS = new File (fileOfData[9], "/taxaName/source/282_panel.txt").getAbsolutePath();
        String panelAmesFileS = new File (fileOfData[9], "/taxaName/source/Ames_January.txt").getAbsolutePath();
        String out282fileS = new File (fileOfData[9], "/taxaName/282.txt").getAbsolutePath();
        String outAmesfileS = new File (fileOfData[9], "/taxaName/Ames.txt").getAbsolutePath();
        
        String taxaNameFileS = new File (fileOfData[9], "/taxaName/aTaxa.txt").getAbsolutePath();
        String out282AmesfileS = new File (fileOfData[9], "/taxaName/282Ames.txt").getAbsolutePath();
        
        String b73FileS = new File (fileOfData[9], "/taxaName/B73.txt").getAbsolutePath();
        
        TaxaNameUtils tnu = new TaxaNameUtils(taxaNameFileS);
        //tnu.mk282List(panel282FileS, out282fileS);
        tnu.mkAmesList(panelAmesFileS, outAmesfileS);
        //tnu.mk282AmesList(panel282FileS, panelAmesFileS, out282AmesfileS);
        //tnu.mkB73InitialList(b73FileS);
    }
    
    public void mergeBinsInBBT () {
        //String binByTaxaNAMFileS = new File (fileOfData[9], "/bbt/cnv.NAM12.mtaxa.bbt.txt").getAbsolutePath();
        //String binByTaxaMBinFileS = new File (fileOfData[9], "/bbt/cnv.NAM12.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        //String binByTaxaFileS = new File (fileOfData[9], "/bbt/cnv.mtaxa.bbt.txt").getAbsolutePath();
        //String binByTaxaMBinFileS = new File (fileOfData[9], "/bbt/cnv.mtaxa.200k.bbt.txt").getAbsolutePath();
        
        String binByTaxaFileS = new File (fileOfData[9], "/bbt/cnv.ataxa.50k.bbt.txt").getAbsolutePath();
        String binByTaxaMBinFileS = new File (fileOfData[9], "/bbt/cnv.ataxa.1000k.bbt.txt").getAbsolutePath();
        
        
        BinByTaxa bbt = new BinByTaxa(binByTaxaFileS);
        bbt.mergeBins(1000000);
        bbt.writeTxtFile(binByTaxaMBinFileS);
    }
    
    public void addB73TaxaInBBT () {
        String binByTaxaFileS = new File (fileOfData[9], "/bbt/cnv.bbt.50k.txt").getAbsolutePath();
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.ataxa.50k.bbt.txt").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaFileS);
        bbt.addb73Taxa();
        //bbt.writeTxtFile(binByTaxaMTaxaFileS);
    }
    
    public void creatTaxaCoverage () {
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/paper.mtaxa.50k.bbt.txt").getAbsolutePath();
        String taxaCoverageFileS = new File (fileOfData[9], "/totalCount/paper.mtaxa.taxaCount.txt").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaMTaxaFileS);
        bbt.mkTaxaCountFile(taxaCoverageFileS);
    }
    
    public void mergedTaxaInBBT () {
        String binByTaxaFileS = new File (fileOfData[9], "/bbt/cnv.bbt.50k.txt").getAbsolutePath();
        String binByTaxaMTaxaFileS = new File (fileOfData[9], "/bbt/cnv.mtaxa.bbt.50k.txt").getAbsolutePath();
        BinByTaxa bbt = new BinByTaxa(binByTaxaFileS);
        bbt.mergeTaxa();
        bbt.writeTxtFile(binByTaxaMTaxaFileS);
    }
    
    public void creatBinByTaxa () {
        String genomeInfoFileS = "E:/Database/InfoFile/ChrLenCentPosi.txt";
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String mergedTBTFileS = "J:/tbt_20.byte";
        BinByTaxa bbt = new BinByTaxa (genomeInfoFileS, tagPMapFileS, mergedTBTFileS);
        String binByTaxaFileS = new File (fileOfData[9], "/bbt/cnv.bbt.50k.bin").getAbsolutePath();
        bbt.writeBinaryFile(binByTaxaFileS);
    }
    
    public void getCNVTagCount () {
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String cnvTagCountFileS =  new File (fileOfData[9], "/tagCount/cnv.cnt").getAbsolutePath();
        PAVUtils util = new PAVUtils();
        util.makeCNVTagCountFile(tagPMapFileS, tagCountFileS, cnvTagCountFileS);
    }
    
    public void getB73TagPosCount () {
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String B73TagPosCountFileS = "E:/Research/pav/attDis/B73TagPosCount.txt";
        String cnvTBTFileS = new File (fileOfData[9], "/tbt/cnv_NAMParents.tbt.byte").getAbsolutePath();
        PAVUtils util = new PAVUtils();
        util.getB73TagPosCount(tagPMapFileS, cnvTBTFileS, B73TagPosCountFileS);
    }
    
    public void identifyPAV () {
        String samFileS = new File (fileOfData[6], "NAM.jpgmap.fil.very-sensitive-local-k10.sam").getAbsolutePath();
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.fil.txt").getAbsolutePath();
        String pavPGMapFileS = new File (fileOfData[8], "pav.jpgmap.fil.txt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap (tagPGMapFileS, false);
        boolean[] ifPAV = pgmap.getIfPAV(samFileS);
        pgmap.writeTxtFile(pavPGMapFileS, ifPAV);
        TagJointPGMap pgmap2 = new TagJointPGMap (pavPGMapFileS, false);
        pgmap2.sortByGeneticPosition();
        pgmap2.writeTxtFile("M:/test.txt");
    }
    
    public void alignTagOfJointPGMap () {
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.fil.txt").getAbsolutePath();
        String fastaFileS = new File (fileOfData[8], "NAM.jpgmap.fil.fas").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap (tagPGMapFileS, false);
        pgmap.writeFastaFile(fastaFileS);
        //align those tags using bowtie2
    }
    
    public void jointStatistics () {
        String filteredTagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.fil.txt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap(filteredTagPGMapFileS, false);
        boolean[] ifB73 = pgmap.getIfRefUniqueTag();
        boolean[] ifNonB73 = pgmap.getReverseFilter(ifB73);
        boolean[] ifSingle = pgmap.getIfSingleMapping(0.0001);
        boolean[] ifMultiple = pgmap.getIfMultipleMapping(0.0001);
        pgmap.getFilterCross(ifB73, ifMultiple);
    }
    
    public void filteredTagJointPGMapByThresh () {
        double pThresh = 0.0001;
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String filteredTagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.fil.txt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap(tagPGMapFileS, false);
        boolean[] ifUnderThresh = pgmap.getIfUnderThresholdEitherFamilyGroup(pThresh);
        pgmap.writeTxtFile(filteredTagPGMapFileS, ifUnderThresh);
    }
    
    public void creatPAVUltraresPosFile () {
        String pavPGMapFileS = new File (fileOfData[8], "pav.fpgmap.ultrares.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String pavPosFileS = "E:/Research/pav/attDis/pavPos.ultrares.txt";
        TagFPGMap fpgmap = new TagFPGMap (pavPGMapFileS);
        fpgmap.writePAVPosFileS(pavPosFileS, tagCountFileS);
    }
    
    public void creatPAVPosFile () {
        String pavPGMapFileS = new File (fileOfData[8], "pav.fpgmap.res.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        String pavPosFileS = "E:/Research/pav/attDis/pavPos.res.txt";
        TagFPGMap fpgmap = new TagFPGMap (pavPGMapFileS);
        fpgmap.writePAVPosFileS(pavPosFileS, tagCountFileS);
    }
    
    public void creatUltraResPAVFPGMap () {
        String finalPGMapFileS = new File (fileOfData[8], "All.fpgmap.ultrares.txt").getAbsolutePath();
        String fastaFileS = new File (fileOfData[13], "fpg.ultrares.fasta").getAbsolutePath();
        TagFPGMap fpgmap = new TagFPGMap (finalPGMapFileS);
        //fpgmap.writeFastaFile(fastaFileS);
        
        String samFileS = new File (fileOfData[6], "fpg.ultrares.very-sensitive-local-k5.sam").getAbsolutePath();
        String pavPGMapFileS = new File (fileOfData[8], "pav.fpgmap.ultrares.txt").getAbsolutePath();
        boolean[] ifPAV = fpgmap.getIfPAV(samFileS);
        fpgmap.writeTxtFile(pavPGMapFileS, ifPAV);
    }
    
    public void sortFPGPaparFile () {
        String pavFPGPaperFileS = "M:/pav/PhyGenMapping/v1.paper.txt";
        String pavFPGPaperSortFileS = "M:/pav/PhyGenMapping/v1.sort.paper.txt";
        PAVUtils util = new PAVUtils();
        util.sortPavFPGPaperFile(pavFPGPaperFileS, pavFPGPaperSortFileS);
    }
    
    
    public void mkPredictionValueFile () {
        String predictionGFileS = new File (fileOfData[12], "g.pre.txt").getAbsolutePath();
        String predictionGJFileS = new File (fileOfData[12], "gj.pre.txt").getAbsolutePath();
        String predictionJFileS = new File (fileOfData[12], "j.pre.txt").getAbsolutePath();
        
        String tagGPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
        String tagGPGMapResFileS = new File (fileOfData[8], "All.gpgmap.res.txt").getAbsolutePath();
        
        String tagGJPGMapFileS = new File (fileOfData[8], "All.mpgmap.txt").getAbsolutePath();
        String tagGJPGMapResFileS = new File (fileOfData[8], "All.mpgmap.res.txt").getAbsolutePath();
        
        String tagJPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String tagJPGMapResFileS = new File (fileOfData[8], "NAM.jpgmap.res.txt").getAbsolutePath();
        
        String gpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.gpgmap.res.predictionValue.txt";
        String gjpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.mpgmap.res.predictionValue.txt";
        String jpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.jpgmap.res.predictionValue.txt";
        
        PAVUtils util = new PAVUtils();
        util.mkGPredictionValueFile(predictionGFileS, tagGPGMapFileS, tagGPGMapResFileS, gpgMapResPredictionValueFileS);
        util.mkGJPredictionValueFile(predictionGJFileS, tagGJPGMapFileS, tagGJPGMapResFileS, gjpgMapResPredictionValueFileS);
        util.mkJPredictionValueFile(predictionJFileS, tagJPGMapFileS, tagJPGMapResFileS, jpgMapResPredictionValueFileS);
    }
    
    public void creatPAVFPGMap () {
        String finalPGMapFileS = new File (fileOfData[8], "All.fpgmap.looseRes.txt").getAbsolutePath();
        String fastaFileS = new File (fileOfData[13], "fpg.fasta").getAbsolutePath();
        TagFPGMap fpgmap = new TagFPGMap (finalPGMapFileS);
        //fpgmap.writeFastaFile(fastaFileS); //bowtie2 align
        
        String TOGMFileS = new File (fileOfData[8], "v1.loose.togm.txt").getAbsolutePath();
        fpgmap.writeTextTOGMFile(TOGMFileS);
        
        String samFileS = new File (fileOfData[6], "fpg.very-sensitive-local-k10.sam").getAbsolutePath();
        String alignmentCheckFileS = new File (fileOfData[8], "fpg.alignCheck.txt").getAbsolutePath();
        String alignmentCheck2FileS = new File (fileOfData[8], "fpg.alignCheck2.txt").getAbsolutePath();
        String pavPGMapFileS = new File (fileOfData[8], "pav.fpgmap.res.txt").getAbsolutePath();
        //fpgmap.writeAlignmentCheck(samFileS, 10000000, alignmentCheckFileS); // upgrade to check2
        //fpgmap.writeAlignmentCheck2(samFileS, 10000000, alignmentCheck2FileS);
        
        //boolean[] ifPAV = fpgmap.getIfPAV(samFileS);
        //fpgmap.writeTxtFile(pavPGMapFileS, ifPAV);
        
        //fpgmap.screenPrintProportionOfMappedTags(samFileS);
        
        //String pavFPGPaperFileS = new File (fileOfData[8], "v1.paper.txt").getAbsolutePath();
        //fpgmap.writeFPGPaperFile(pavFPGPaperFileS, ifPAV);
    }
    
    public void mergeHighResPGMap () {
        String highResGJFileS = new File (fileOfData[8], "All.mpgmap.res.txt").getAbsolutePath();
        String highResGFileS = new File (fileOfData[8], "All.gpgmap.res.txt").getAbsolutePath();
        String highResJFileS = new File (fileOfData[8], "NAM.jpgmap.res.txt").getAbsolutePath();
        String finalPGFileS = new File (fileOfData[8], "All.fpgmap.res.txt").getAbsolutePath();
        
        String gpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.gpgmap.res.predictionValue.txt";
        String gjpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.mpgmap.res.predictionValue.txt";
        String jpgMapResPredictionValueFileS = "M:/pav/PhyGenMapping/predictionValue/All.jpgmap.res.predictionValue.txt";
        //TagFPGMap fpgmap = new TagFPGMap (highResGJFileS, highResGFileS, highResJFileS, gjpgMapResPredictionValueFileS, gpgMapResPredictionValueFileS, jpgMapResPredictionValueFileS);
        //fpgmap.checkMappingQuality();
        //fpgmap.writeTxtFile(finalPGFileS);
        
        TagFPGMap fpgmap = new TagFPGMap (finalPGFileS);
        fpgmap.checkMappingQuality();
        /**********Deprecated***********/
        //String highResGJFileS = new File (fileOfData[8], "All.mpgmap.ultrares.txt").getAbsolutePath();
        //String highResGFileS = new File (fileOfData[8], "All.gpgmap.ultrares.txt").getAbsolutePath();
        //String finalPGFileS = new File (fileOfData[8], "All.fpgmap.ultrares.txt").getAbsolutePath();
        //TagFPGMap fpgmap = new TagFPGMap (highResGJFileS, highResGFileS);
        //fpgmap.checkMappingQuality();
        //fpgmap.writeTxtFile(finalPGFileS);
        /*******************************/
    }
    
    public void updateResolution (int type) {
        if (type == 1) {
            String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "g.pre.txt").getAbsolutePath();
            String highResFileS = new File (fileOfData[8], "All.gpgmap.res.txt").getAbsolutePath();
            TagGWASPGMap gpgmap = new TagGWASPGMap(tagGWASPGMapFileS);
            //boolean[] ifOut = gpgmap.getIfHighResolution(predictionFileS, Math.log10(50000));
            //gpgmap.writeTxtFile(highResFileS, ifOut);
            
            //String ultraHighResFileS = new File (fileOfData[8], "All.gpgmap.ultrares.txt").getAbsolutePath();
            //boolean[] ifOut = gpgmap.getIfHighResolution(predictionFileS, Math.log10(2000));
            //gpgmap.writeTxtFile(ultraHighResFileS, ifOut);
            
            String looseResFileS = new File (fileOfData[8], "All.gpgmap.looseRes.txt").getAbsolutePath();
            boolean[] ifOut = gpgmap.getIfHighResolution(predictionFileS, Math.log10(2000000));
            gpgmap.writeTxtFile(looseResFileS, ifOut);
        }
        else if (type == 2) {
            String mergedPGMapFileS = new File (fileOfData[8], "All.mpgmap.txt").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "gj.pre.txt").getAbsolutePath();
            String highResFileS = new File (fileOfData[8], "All.mpgmap.res.txt").getAbsolutePath();
            TagMPGMap mpgmap = new TagMPGMap (mergedPGMapFileS);
            //boolean[] ifOut = mpgmap.getIfHighResolution(predictionFileS, Math.log10(100000));
            //mpgmap.writeTxtFile(highResFileS, ifOut);
            
            //String ultraHighResFileS = new File (fileOfData[8], "All.mpgmap.ultrares.txt").getAbsolutePath();
            //boolean[] ifOut = mpgmap.getIfHighResolution(predictionFileS, Math.log10(2000));
            //mpgmap.writeTxtFile(ultraHighResFileS, ifOut);
            
            String looseResFileS = new File (fileOfData[8], "All.mpgmap.looseRes.txt").getAbsolutePath();
            boolean[] ifOut = mpgmap.getIfHighResolution(predictionFileS, Math.log10(10000000));
            mpgmap.writeTxtFile(looseResFileS, ifOut);
        }
        else if (type == 3) {
            String tagJPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "j.pre.txt").getAbsolutePath();
            String highResFileS = new File (fileOfData[8], "NAM.jpgmap.res.txt").getAbsolutePath();
            TagJointPGMap jpgmap = new TagJointPGMap(tagJPGMapFileS, false);
            boolean[] ifOut = jpgmap.getIfHighResolution(predictionFileS, Math.log10(100000));
            jpgmap.writeTxtFile(highResFileS, ifOut);
            
            //String looseResFileS = new File (fileOfData[8], "NAM.jpgmap.looseRes.txt").getAbsolutePath();
            //boolean[] ifOut = jpgmap.getIfHighResolution(predictionFileS, Math.log10(5000000));
            //jpgmap.writeTxtFile(looseResFileS, ifOut);
        }
    }
    
    public void creatPredictFile (int type) {
        String recombinationFileS = new File (fileOfData[10], "recombination.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        if (type == 1) {
            String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
            String wekaFileS = new File (fileOfData[12], "g.arff").getAbsolutePath();
            TagGWASPGMap gpgmap = new TagGWASPGMap(tagGWASPGMapFileS);
            gpgmap.mkWekaFile(tagCountFileS, recombinationFileS, wekaFileS);
            
            String modelFileS = new File (fileOfData[10], "train_result/G_model.model").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "g.pre.txt").getAbsolutePath();
            Weka w = new Weka(wekaFileS, modelFileS, predictionFileS);
        }
        else if (type == 2) {
            String mergedPGMapFileS = new File (fileOfData[8], "All.mpgmap.txt").getAbsolutePath();
            //String mergedPGMapFileS = new File (fileOfData[10], "UB73train.mpgmap.txt").getAbsolutePath();
            String wekaFileS = new File (fileOfData[12], "gj.arff").getAbsolutePath();
            TagMPGMap mpgmap = new TagMPGMap (mergedPGMapFileS);
            mpgmap.mkWekaFile(tagCountFileS, recombinationFileS, wekaFileS);
            
            String modelFileS = new File (fileOfData[10], "train_result/GJ_model.model").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "gj.pre.txt").getAbsolutePath();
            Weka w = new Weka(wekaFileS, modelFileS, predictionFileS);
        }
        else if (type == 3) {
            String tagJPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
            String wekaFileS = new File (fileOfData[12], "j.arff").getAbsolutePath();
            TagJointPGMap jpgmap = new TagJointPGMap(tagJPGMapFileS, false);
            jpgmap.mkWekaFile(tagCountFileS, recombinationFileS, wekaFileS);
            
            String modelFileS = new File (fileOfData[10], "train_result/J_model.model").getAbsolutePath();
            String predictionFileS = new File (fileOfData[12], "j.pre.txt").getAbsolutePath();
            Weka w = new Weka(wekaFileS, modelFileS, predictionFileS);
        }
    }
    
    public void mergePGMaps () {
        String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
        String tagJPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String mergedPGMapFileS = new File (fileOfData[8], "All.mpgmap.txt").getAbsolutePath();
        TagMPGMap mpgmap = new TagMPGMap(tagGWASPGMapFileS, tagJPGMapFileS);
        mpgmap.writeTxtFile(mergedPGMapFileS);
    }
    
    public void checkGWASMappingQuality () {
        String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
        TagGWASPGMap gpgmap = new TagGWASPGMap (tagGWASPGMapFileS);
        gpgmap.checkGMappingQuality(0.01, 10);
    }
    
    public void creatDiJointMappingFile () {
        String tagPGMapFileS = new File (fileOfData[14], "NAM.jpgmap.uniqueGene.txt").getAbsolutePath();
        String tagDiPGMapGeneFileS = new File (fileOfData[14], "NAM.jpgmap.di.gene.txt").getAbsolutePath();
        
        TagJointPGMap pgmap = new TagJointPGMap(tagPGMapFileS, false);
        pgmap.sortByPhysicalPosition();
        pgmap.orderFamilyGroupByReference();
        boolean[] ifDi = pgmap.getIfMultipleMapping(0.00001);
        pgmap.writeTxtFile(tagDiPGMapGeneFileS, ifDi);
        
        tagPGMapFileS = new File (fileOfData[14], "NAM.jpgmap.uniqueNonGene.txt").getAbsolutePath();
        String tagDiPGMapNonGeneFileS = new File (fileOfData[14], "NAM.jpgmap.di.nongene.txt").getAbsolutePath();
        pgmap = new TagJointPGMap(tagPGMapFileS, false);
        pgmap.sortByPhysicalPosition();
        pgmap.orderFamilyGroupByReference();
        ifDi = pgmap.getIfMultipleMapping(0.00001);
        pgmap.writeTxtFile(tagDiPGMapNonGeneFileS, ifDi);
    }
    
    
    public void creatGeneUniqueJointMappingFile () {
        String tagUniquePGMapFileS = new File (fileOfData[14], "NAM.jpgmap.unique.txt").getAbsolutePath();
        String fastaFileS = new File (fileOfData[14], "NAM.jpgmap.unique.fasta").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap(tagUniquePGMapFileS, false);
        //pgmap.writeFastaFile(fastaFileS);
        
        String samFileS = new File (fileOfData[14], "jpgmap.unique-gene-k2.sam").getAbsolutePath();
        String tagUniqueGenePGMapFileS = new File (fileOfData[14], "NAM.jpgmap.uniqueGene.txt").getAbsolutePath();
        String tagUniqueNonGenePGMapFileS = new File (fileOfData[14], "NAM.jpgmap.uniqueNonGene.txt").getAbsolutePath();
        boolean[] ifMatchUniqueGene = pgmap.getIfMatchUniqueGene(samFileS) ;
        pgmap.writeTxtFile(tagUniqueGenePGMapFileS, ifMatchUniqueGene);
        
        boolean[] ifMatchUniqueNonGene = pgmap.getIfMatchUniqueNonGene(samFileS) ;
        pgmap.writeTxtFile(tagUniqueNonGenePGMapFileS, ifMatchUniqueNonGene);
    }
    
    public void creatUniqueJointMappingFile () {
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String tagUniquePGMapFileS = new File (fileOfData[14], "NAM.jpgmap.unique.txt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap(tagPGMapFileS, false);
        boolean[] ifUnique = pgmap.getIfRefUniqueTag();
        pgmap.writeTxtFile(tagUniquePGMapFileS, ifUnique);
    }
    
    public void checkJointMappingQuality () {
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap(tagPGMapFileS, false);
        pgmap.checkJMappingQuality(0.0001, true, true);
        //pgmap.checkMappingQualityMultiGroup(0.0001, 1);
    }
    
    public void creatTagGWASPGMap () {
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String tagGWASGMapFileS = new File (fileOfData[5], "merged.mapping.gwas.txt").getAbsolutePath();
        String tagGWASPGMapFileS = new File (fileOfData[8], "All.gpgmap.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        TagGWASPGMap gpgmap = new TagGWASPGMap (tagPMapFileS, tagGWASGMapFileS, tagCountFileS);
        gpgmap.writeTxtFile(tagGWASPGMapFileS);
    }
    
    
    public void creatTagJointPGMap () {
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k2.pmap").getAbsolutePath();
        String tagGMapFileS = new File (fileOfData[5], "merged.mapping.joint.txt").getAbsolutePath();
        String tagPGMapFileS = new File (fileOfData[8], "NAM.jpgmap.txt").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        TagJointPGMap pgmap = new TagJointPGMap (tagPMapFileS, tagGMapFileS, tagCountFileS);
        pgmap.writeTxtFile(tagPGMapFileS);
    }
    
    public void mapTBTJoint () {
        String tbtFileS = new File (fileOfData[1], "tbt_test.byte").getAbsolutePath();
        String anchorMapFileS = new File (fileOfData[4], "testData/HMM/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_HMM_filtSites_imp_chr+.hmp.txt").getAbsolutePath();
        //String anchorMapFileS = new File (fileOfData[4], "impNAM/NAMc+.hmp.txt").getAbsolutePath();
        String pedigreeFileS = new File (fileOfData[7], "NAM.txt").getAbsolutePath();
        String mappingResultFileS = "M:/mapR.txt";
        TagAgainstAnchorNAM taa = new TagAgainstAnchorNAM (tbtFileS, anchorMapFileS, pedigreeFileS, mappingResultFileS);
        
    }
    
    public void mapTBTSorghum () {
        String tbtFileS = "M:/sorghum/test.tbt.byte";
        String anchorMapFileS = "M:/sorghum/chr+.hmp.txt";
        //String anchorMapFileS = new File (fileOfData[4], "impNAM/NAMc+.hmp.txt").getAbsolutePath();
        String pedigreeFileS = "M:/sorghum/SbPedigree.txt";
        String mappingResultFileS = "M:/mapR.txt";
        TagAgainstAnchorBiParental taa = new TagAgainstAnchorBiParental (tbtFileS, anchorMapFileS, pedigreeFileS, mappingResultFileS);
        
    }
    
    public void selectNamFromAllAnchor () {
        String pedigreeFileS = new File (fileOfData[7], "NAM.txt").getAbsolutePath();
        String inAnchorMapDirS = new File (fileOfData[4], "impAll").getAbsolutePath();
        String outAnchorMapDirS = new File (fileOfData[4], "impNAM").getAbsolutePath();
        Pedigree ped = new Pedigree (pedigreeFileS);
        
        
        String[][] samples = ped.getSampleOfFamilies(ped.getFamilyStartWith("NAM_", 100), 100);
        File[] inHapMapFiles = new File (inAnchorMapDirS).listFiles();
        for (int i = 0; i < inHapMapFiles.length; i++) {
            String[] temp = inHapMapFiles[i].getName().split("\\.");
            String tem = "NAM" + temp[temp.length-3];
            String outfileS = outAnchorMapDirS + "/" + tem + ".hmp.txt";
            HapMapUtils hmu = new HapMapUtils(inHapMapFiles[i].getAbsolutePath(), outfileS);
            hmu.checkSegregationInFamilyAndOutput(samples, (float)0.2);
            System.out.println("BiParental families selected from " + inHapMapFiles[i]);
        }
    }
    
    public void creatTagPMap () {
        String samFileS = new File (fileOfData[6], "26Mtag-very-sensitive-local-k20.sam").getAbsolutePath();
        String tagPMapFileS = new File (fileOfData[6], "26Mtag-k20.pmap").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        new TagPMap(tc, samFileS).writeTagPMap(tagPMapFileS);
    }
    
    public void creatTagMap () {
        //String gMappingFileS = new File (fileOfData[5], "merged.mapping.txt").getAbsolutePath();
        String gMappingFileS = "M:/sliceTBT-1000.mapping.txt";
        String samFileS = new File (fileOfData[6], "26Mtag-very-sensitive-local-k3.sam").getAbsolutePath();
        String tagCountFileS = new File (fileOfData[2], "merged_20.cnt").getAbsolutePath();
        //String gMappingUpdateFileS = new File (fileOfData[5], "update.mapping.txt").getAbsolutePath();
        String gMappingUpdateFileS = "M:/updateB73mapR.txt";
        TagMap_Deprecated tm = new TagMap_Deprecated (gMappingFileS, samFileS, tagCountFileS);
        tm.checkLdAndBlast();
        tm.updateMappingFile(gMappingFileS, gMappingUpdateFileS, true);
        
    }
    
    public void mergeMappingResult (boolean ifJointLinkage) {
        if (ifJointLinkage) {
            String inputDirS = new File (fileOfData[5], "joint").getAbsolutePath();
            String mergedFileS = fileOfData[5].getAbsolutePath() + "/merged.mapping.joint.txt";
            new PAVUtils().mergeGeneticMappingFiles(inputDirS, mergedFileS);
        }
        else {
            String inputDirS = new File (fileOfData[5], "gwas").getAbsolutePath();
            String mergedFileS = fileOfData[5].getAbsolutePath() + "/merged.mapping.gwas.txt";
            new PAVUtils().mergeGeneticMappingFiles(inputDirS, mergedFileS);
        }
    }
    
    public void sliceTBT () {
        String sliceTBTDirS = fileOfData[3].getAbsolutePath();
        //String mergedTBT = new File (fileOfData[1], "434GFAAXX_s_2.tbt.byte").getAbsolutePath();
        String mergedTBT = "Z:/tbt_20.byte";
        int sliceNum = 200;
        PAVUtils uti = new PAVUtils();
        uti.sliceTBT(sliceTBTDirS, mergedTBT, sliceNum);
    }
    
    
    public void mapTBTGWAS () {//this step is just a test, mapping runs on clusters,
        String tbtFileS = new File (fileOfData[1], "tbt_test_50.byte").getAbsolutePath();
        //String tbtFileS = "M:/test55againstAnchor/fake1.tbt.byte";
        String anchorMapFileS = new File (fileOfData[4], "testData/impAll/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt").getAbsolutePath();
        //String anchorMapFileS = "N:/Zea/build20120110/imp/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt";
        String mappingResultFileS = "M:/mapR.txt";
        //TagCallerAgainstAnchorMT tcaa=new TagCallerAgainstAnchorMT(tbtFileS, anchorMapFileS, null, mappingResultFileS, -1, -1);
        TagAgainstAnchor taa = new TagAgainstAnchor (tbtFileS, anchorMapFileS, mappingResultFileS);
    }
    
    public void mergeTBTByRow () {
        MergeTagsByTaxaFilesByRowPlugin mtbt = new MergeTagsByTaxaFilesByRowPlugin ();
        mtbt.creatMergeTBTByTagCount(fileOfData[0].getAbsolutePath(), new File (fileOfData[1], "tbt_20.byte").getAbsolutePath(), new File (fileOfData[2], "merged_20.cnt").getAbsolutePath());
    }
    
    public void getTagPairs () {
        String mergedTagCountOfAllS = fileOfData[2].getAbsolutePath() + "/merged_20.cnt";
        String tagPairS = fileOfData[3].getAbsolutePath() + "/tagPair.tps";
        double etr = 0.05;
        UNetworkFilter unf = new UNetworkFilter (mergedTagCountOfAllS, etr, tagPairS);
    }
    
    public void creatDir (String workingDirS) {
        fileOfData = new File[childDir.length];
        for (int i = 0; i < childDir.length; i++) {
            fileOfData[i] = new File(parentDir, childDir[i]);
            fileOfData[i].mkdirs();
        }
    }
    
    public static void main (String[] args) {
        String workingDirS = "M:/pav/";
        new PAVgo (workingDirS);
    }
}
