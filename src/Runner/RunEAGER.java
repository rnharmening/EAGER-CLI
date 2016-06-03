/*
 * Copyright (c) 2016. EAGER-CLI Alexander Peltzer
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package Runner;

import IO.Communicator;
import IO.FileSearcher;
import Modules.filehandling.*;
import Modules.genotyping.*;
import Modules.indexing.*;
import Modules.mapping.*;
import Modules.preprocessing.ClipAndMerge;
import Modules.preprocessing.FastQC;
import Modules.preprocessing.QualityTrimmer;
import Modules.stats.*;
import com.thoughtworks.xstream.XStream;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by peltzer on 22.01.14.
 */
public class RunEAGER {
    private Communicator communicator;
    private ArrayList<ModulePool> pools = new ArrayList<ModulePool>();


    public RunEAGER(Communicator communicator) throws IOException, InterruptedException {
        this.communicator = communicator;
        checkForConfigured();
        executeAll();
    }

    public static void main(String[] args) throws IOException, ClassNotFoundException, InterruptedException {
        if (args.length < 1) {
            System.err.println("Please specify [n] input files to be used in the EAGER Pipeline. This could be for example: \n");
            System.err.println(" - a folder with subfolders containing EAGER Configuration files, which will be then searched and executed.");
            System.err.println(" - a single EAGER Configuration file, created using the GUI.");
            System.err.println(" - (optional) a properties file, describing paths to tools required by the pipeline.");
            System.exit(0);
        } else {
            ArrayList<String> inputFiles = new ArrayList<String>();
            FileSearcher fileSearcher = new FileSearcher();
            inputFiles = fileSearcher.processFiles(args[0]);
            System.out.println("Found " + inputFiles.size() + " input configuration files.");
            int counter = 0;

            for (String inputFile : inputFiles) {
                counter++;
                System.out.println("Processing file # " + counter);
                System.out.println("Schaffa, Schaffa, Genome baua!");
                XStream xstream = new XStream();
                InputStream in = new FileInputStream(inputFile);
                Communicator c = (Communicator) xstream.fromXML(in);
                RunEAGER runEAGER = new RunEAGER(c);
            }

        }
    }

    /**
     * Method checking for proper configuration of the pipeline and which subpipeline should be executed by EAGER.
     *
     * @throws IOException
     * @throws InterruptedException
     */

    private void checkForConfigured() throws IOException, InterruptedException {
        if (this.communicator.isOrganism()) {
            if (this.communicator.isOrganismage()) {
                createAncientHumanPipeline();
            } else {
                createHumanPipeline();
            }
        } else {
            if (this.communicator.isOrganismage()) {
                createAncientBacterialPipeline();
            } else {
                createBacterialPipeline();
            }
        }
    }


    private void executeAll() throws IOException, InterruptedException {
        //Set Input Path for first pool correctly and automatically
        pools.get(0).setCurrentFilePath(communicator.getGUI_inputfiles());
        //Now start all the pools
        for (ModulePool pool : pools) {
            pool.start();
        }
    }


    /**
     * This code here generates the corresponding pipeline Modules depending on the configuration of the
     * Communicator Object.
     */

    /**
     * Creates a Pipeline for modern bacterial input data.
     */

    private void createBacterialPipeline() throws IOException, InterruptedException {
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        ModulePool bacterialpool = new ModulePool();
        ModulePool preprocesspool = new ModulePool();
        ModulePool gatkpool = new ModulePool();
        preprocesspool.addModule(new CreateResultsDirectories(communicator));
        if (communicator.isReferenceselected()) {
            preprocesspool.addModule(new ReferenceRenamer(communicator));
        }

        if (communicator.isRun_fastqc()) {
            preprocesspool.addModule(new FastQC(communicator));
        }

        if (communicator.isRun_clipandmerge()) {
            if(communicator.getMerge_type().equals("PAIRED")){
                if (!communicator.isMerge_only_clipping()) {
                    preprocesspool.addModule(new ClipAndMerge(communicator));
                } else {
                    preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.ADAPTER_CLIPPING_ONLY));
                }
            } else {
               preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.SINGLE_ENDED_ONLY));
            }
        }

        if (communicator.isRun_qualityfilter()) {
            preprocesspool.addModule(new QualityTrimmer(communicator));
            preprocesspool.addModule(new FastQC(communicator, FastQC.AFTERMERGING));
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("CircularMapper")) {
            bacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                bacterialpool.addModule(new BWAIndex(communicator));
            }
            addCircularMapping(bacterialpool);
        }


        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWA")) {
            bacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                bacterialpool.addModule(new BWAIndex(communicator));
            }
            addBWAMapping(bacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Bowtie 2")) {
            bacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("BT2")) {
                bacterialpool.addModule(new BT2Index(communicator));
            }
            addBT2Mapping(bacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Stampy")) {
            bacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("Stampy")) {
                bacterialpool.addModule(new StampyIndex(communicator));
                bacterialpool.addModule(new StampyHash(communicator));
            }
            addStampyMapping(bacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWAMem")) {
            bacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                bacterialpool.addModule(new BWAIndex(communicator));
            }
            addBWAMemMapping(bacterialpool);
        }

        if (communicator.isRun_complexityestimation()) {
            addComplexityEstimation(bacterialpool);
        }

        if (communicator.isRmdup_run() && !communicator.isMarkdup_run()) {
            bacterialpool.addModule(new DeDup(communicator));
        }

        if (communicator.isMarkdup_run()) {
            bacterialpool.addModule(new MarkDuplicates(communicator));
        }

        if(communicator.isRun_mapping() && (communicator.isMarkdup_run() || communicator.isRmdup_run())){
            bacterialpool.addModule(new SamtoolsIndex(communicator, SamtoolsIndex.DEDUP));
        }

        if(communicator.isSchmutzi_run()){
            addContaminationEstimation(bacterialpool);
        }


        if (communicator.isRun_coveragecalc() && communicator.isRun_mt_capture_mode()) {
            bacterialpool.addModule(new QualiMap(communicator,QualiMap.CAPTURE));
        } else if (communicator.isRun_coveragecalc()) {

            bacterialpool.addModule(new QualiMap(communicator));
        }

        if (communicator.isRun_mapdamage()) {
            bacterialpool.addModule(new MapDamage(communicator));
        }

        if (communicator.isRun_gatksnpcalling() && !communicator.getGatk_caller().equals("ANGSD")) {
            gatkpool = createGATKSNPCallingPipeline();
            gatkpool.setCurrentFilePath(bacterialpool.getCurrentFilePath());
        }

        if(communicator.isRun_gatksnpcalling() && communicator.getGatk_caller().equals("ANGSD")){
            gatkpool = createANGSDCallingPipeline();
            gatkpool.setCurrentFilePath(bacterialpool.getCurrentFilePath());
        }

        if (communicator.isRun_gatksnpfiltering()) {
            gatkpool.addModule(new GATKVariantFilter(communicator));
            if (!communicator.isRun_vcf2draft()) {
                gatkpool.addModule(new BGZip(communicator));
                gatkpool.addModule(new Tabix(communicator));
            }
        }

        if (communicator.isRun_vcf2draft()) {
            gatkpool.addModule(new VCF2Genome(communicator));
        }

        if (communicator.isRun_cleanup()) {
            gatkpool.addModule(new CleanUpRedundantData(communicator));

        }

        if (communicator.isRun_reportgenerator()) {
            gatkpool.addModule(new ReportGenerator(communicator));
        }

        pools.add(preprocesspool);
        bacterialpool.addPredecessor(preprocesspool);
        pools.add(bacterialpool);
        gatkpool.addPredecessor(bacterialpool);
        pools.add(gatkpool);

    }

    /**
     * Creates a pipeline for ancient bacterial input data
     */

    private void createAncientBacterialPipeline() {
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        ModulePool preprocesspool = new ModulePool();
        ModulePool ancientbacterialpool = new ModulePool();
        ModulePool gatkpool = new ModulePool();

        preprocesspool.addModule(new CreateResultsDirectories(communicator));
        if (communicator.isReferenceselected()) {
            preprocesspool.addModule(new ReferenceRenamer(communicator));
        }

        if (communicator.isRun_fastqc()) {
            preprocesspool.addModule(new FastQC(communicator));
        }
        if (communicator.isRun_clipandmerge()) {
            if(communicator.getMerge_type().equals("PAIRED")){
                if (!communicator.isMerge_only_clipping()) {
                    preprocesspool.addModule(new ClipAndMerge(communicator));
                } else {
                    preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.ADAPTER_CLIPPING_ONLY));
                }
            } else {
                preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.SINGLE_ENDED_ONLY));
            }
        }
        if (communicator.isRun_qualityfilter()) {
            preprocesspool.addModule(new QualityTrimmer(communicator));
            preprocesspool.addModule(new FastQC(communicator, FastQC.AFTERMERGING));
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("CircularMapper")) {
            ancientbacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                ancientbacterialpool.addModule(new BWAIndex(communicator));
            }
            addCircularMapping(ancientbacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWA")) {
            ancientbacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                ancientbacterialpool.addModule(new BWAIndex(communicator));
            }
            addBWAMapping(ancientbacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Stampy")) {
            ancientbacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("Stampy")) {
                ancientbacterialpool.addModule(new StampyIndex(communicator));
                ancientbacterialpool.addModule(new StampyHash(communicator));
            }
            addStampyMapping(ancientbacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Bowtie 2")) {
            ancientbacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("BT2")) {
                ancientbacterialpool.addModule(new BT2Index(communicator));
            }
            addBT2Mapping(ancientbacterialpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWAMem")) {
            ancientbacterialpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                ancientbacterialpool.addModule(new BWAIndex(communicator));
            }
            addBWAMemMapping(ancientbacterialpool);
        }

        if (communicator.isRun_complexityestimation()) {
            ancientbacterialpool.addModule(new PreseqCCurveCalculation(communicator));
            ancientbacterialpool.addModule(new PreseqLCExtrapCalculation(communicator));
        }

        if (communicator.isRmdup_run() && !communicator.isMarkdup_run()) {
            ancientbacterialpool.addModule(new DeDup(communicator));
        }

        if (communicator.isMarkdup_run()) {
            ancientbacterialpool.addModule(new MarkDuplicates(communicator));
        }

        if(communicator.isRun_mapping() && (communicator.isMarkdup_run() || communicator.isRmdup_run())){
            ancientbacterialpool.addModule(new SamtoolsIndex(communicator, SamtoolsIndex.DEDUP));
        }

        if(communicator.isSchmutzi_run()){
            addContaminationEstimation(ancientbacterialpool);
        }

        if (communicator.isRun_coveragecalc() && communicator.isRun_mt_capture_mode()) {
            ancientbacterialpool.addModule(new QualiMap(communicator,QualiMap.CAPTURE));
        } else if (communicator.isRun_coveragecalc()) {
            ancientbacterialpool.addModule(new QualiMap(communicator));
        }

        if (communicator.isRun_mapdamage()) {
            ancientbacterialpool.addModule(new MapDamage(communicator));
            //We dont run Mapdamage on CM results!
        }

        if (communicator.isRun_gatksnpcalling() && !communicator.getGatk_caller().equals("ANGSD")) {
            gatkpool = createGATKSNPCallingPipeline();
            gatkpool.setCurrentFilePath(ancientbacterialpool.getCurrentFilePath());
        }

        if(communicator.isRun_gatksnpcalling() && communicator.getGatk_caller().equals("ANGSD")){
            gatkpool = createANGSDCallingPipeline();
            gatkpool.setCurrentFilePath(ancientbacterialpool.getCurrentFilePath());
        }

        if (communicator.isRun_gatksnpfiltering()) {
            gatkpool.addModule(new GATKVariantFilter(communicator));
            if (!communicator.isRun_vcf2draft()) {
                gatkpool.addModule(new BGZip(communicator));
                gatkpool.addModule(new Tabix(communicator));
            }
        }

        if (communicator.isRun_vcf2draft()) {
            gatkpool.addModule(new VCF2Genome(communicator));
        }

        if (communicator.isRun_cleanup()) {
            gatkpool.addModule(new CleanUpRedundantData(communicator));
        }

        if (communicator.isRun_reportgenerator()) {
            gatkpool.addModule(new ReportGenerator(communicator));
        }

        pools.add(preprocesspool);
        ancientbacterialpool.addPredecessor(preprocesspool);
        pools.add(ancientbacterialpool);
        gatkpool.addPredecessor(ancientbacterialpool);
        pools.add(gatkpool);

    }

    /**
     * Creates a pipeline for modern human data
     */

    private void createHumanPipeline() {
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        ModulePool preprocesspool = new ModulePool();
        ModulePool humanmodernpool = new ModulePool();
        ModulePool gatkpool = new ModulePool();
        ModulePool reportpool = new ModulePool();

        preprocesspool.addModule(new CreateResultsDirectories(communicator));
        if (communicator.isReferenceselected()) {
            preprocesspool.addModule(new ReferenceRenamer(communicator));
        }


        if (communicator.isRun_fastqc()) {
            preprocesspool.addModule(new FastQC(communicator));
        }
        if (communicator.isRun_clipandmerge()) {
            if(communicator.getMerge_type().equals("PAIRED")){
                if (!communicator.isMerge_only_clipping()) {
                    preprocesspool.addModule(new ClipAndMerge(communicator));
                } else {
                    preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.ADAPTER_CLIPPING_ONLY));
                }
            } else {
                preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.SINGLE_ENDED_ONLY));
            }
        }
        if (communicator.isRun_qualityfilter()) {
            preprocesspool.addModule(new QualityTrimmer(communicator));
            preprocesspool.addModule(new FastQC(communicator, FastQC.AFTERMERGING));
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWA")) {
            humanmodernpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanmodernpool.addModule(new BWAIndex(communicator));
            }
            addBWAMapping(humanmodernpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Stampy")) {
            humanmodernpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("Stampy")) {
                humanmodernpool.addModule(new StampyIndex(communicator));
                humanmodernpool.addModule(new StampyHash(communicator));
            }
            addStampyMapping(humanmodernpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWAMem")) {
            humanmodernpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanmodernpool.addModule(new BWAIndex(communicator));
            }
            addBWAMemMapping(humanmodernpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Bowtie 2")) {
            humanmodernpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("BT2")) {
                humanmodernpool.addModule(new BT2Index(communicator));
            }
            addBT2Mapping(humanmodernpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("CircularMapper")) {
            humanmodernpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanmodernpool.addModule(new BWAIndex(communicator));
            }
            addCircularMapping(humanmodernpool);
        }

        if (communicator.isRun_complexityestimation()) {
            addComplexityEstimation(humanmodernpool);
        }

        if (communicator.isRmdup_run() && !communicator.isMarkdup_run()) {
            humanmodernpool.addModule(new DeDup(communicator));
        }

        if (communicator.isMarkdup_run()) {
            humanmodernpool.addModule(new CleanSam(communicator));
            humanmodernpool.addModule(new MarkDuplicates(communicator));
        }

        if(communicator.isRun_mapping() && (communicator.isMarkdup_run() || communicator.isRmdup_run())){
            humanmodernpool.addModule(new SamtoolsIndex(communicator, SamtoolsIndex.DEDUP));
        }

        if(communicator.isSchmutzi_run()){
            addContaminationEstimation(humanmodernpool);
        }


        if (communicator.isRun_coveragecalc()) {
            if(communicator.isSnpcapturedata() || communicator.isRun_mt_capture_mode()){
                humanmodernpool.addModule(new QualiMap(communicator,QualiMap.CAPTURE));
            } else {
                humanmodernpool.addModule(new QualiMap(communicator));
            }
            humanmodernpool.addModule(new MTToNucRatioCalculator(communicator));
        }


        if (communicator.isRun_gatksnpcalling() && !communicator.getGatk_caller().equals("ANGSD")) {
            gatkpool = createGATKSNPCallingPipeline();
            gatkpool.addPredecessor(humanmodernpool);
            gatkpool.setCurrentFilePath(humanmodernpool.getCurrentFilePath());
        }

        if(communicator.isRun_gatksnpcalling() && communicator.getGatk_caller().equals("ANGSD")){
            gatkpool = createANGSDCallingPipeline();
            gatkpool.addPredecessor(humanmodernpool);
            gatkpool.setCurrentFilePath(humanmodernpool.getCurrentFilePath());
        }

        if (communicator.isRun_gatksnpfiltering()) {
            gatkpool.addModule(new GATKVariantFilter(communicator));
            if (!communicator.isRun_vcf2draft()) {
                gatkpool.addModule(new BGZip(communicator));
                gatkpool.addModule(new Tabix(communicator));
            }
        }

        if (communicator.isRun_vcf2draft()) {
            gatkpool.addModule(new VCF2Genome(communicator));

        }

        if (communicator.isRun_cleanup()) {
            gatkpool.addModule(new CleanUpRedundantData(communicator));
        }

        if (communicator.isRun_reportgenerator()) {
            reportpool.addModule(new ReportGenerator(communicator));
        }

        pools.add(preprocesspool);
        if ( humanmodernpool.getModules().size() > 0 ) {
            humanmodernpool.addPredecessor(preprocesspool);
            pools.add(humanmodernpool);
            reportpool.addPredecessor (humanmodernpool);
        } else {
            System.out.println("Skipping Human Modern Modules");
        }
        if ( gatkpool.getModules().size() > 0 ) {
            gatkpool.addPredecessor(humanmodernpool);
            pools.add(gatkpool);
            reportpool.addPredecessor(gatkpool);
        } else {
            System.out.println("Skipping GATK Modules");
        }
        if ( reportpool.getModules().size() > 0 ) {
            pools.add (reportpool);
        }
    }


    /**
     * Creates a pipeline for ancient human data
     */

    private void createAncientHumanPipeline() throws IOException, InterruptedException {
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        ModulePool preprocesspool = new ModulePool();
        ModulePool humanancientpool = new ModulePool();
        ModulePool gatkpool = new ModulePool();
        ModulePool reportpool = new ModulePool();

        preprocesspool.addModule(new CreateResultsDirectories(communicator));
        if (communicator.isReferenceselected()) {
            preprocesspool.addModule(new ReferenceRenamer(communicator));
        }



        if (communicator.isRun_fastqc()) {
            preprocesspool.addModule(new FastQC(communicator));
        }
        if (communicator.isRun_clipandmerge()) {
            if(communicator.getMerge_type().equals("PAIRED")){
                if (!communicator.isMerge_only_clipping()) {
                    preprocesspool.addModule(new ClipAndMerge(communicator));
                } else {
                    preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.ADAPTER_CLIPPING_ONLY));
                }
            } else {
                preprocesspool.addModule(new ClipAndMerge(communicator, ClipAndMerge.SINGLE_ENDED_ONLY));
            }
        }
        if (communicator.isRun_qualityfilter()) {
            preprocesspool.addModule(new QualityTrimmer(communicator));
            preprocesspool.addModule(new FastQC(communicator, FastQC.AFTERMERGING));

        }
        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWA")) {
            humanancientpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanancientpool.addModule(new BWAIndex(communicator));
            }
            addBWAMapping(humanancientpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Stampy")) {
            humanancientpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("Stampy")) {
                humanancientpool.addModule(new StampyIndex(communicator));
                humanancientpool.addModule(new StampyHash(communicator));
            }
            addStampyMapping(humanancientpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("Bowtie 2")) {
            humanancientpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("BT2")) {
                humanancientpool.addModule(new BT2Index(communicator));
            }
            addBT2Mapping(humanancientpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("BWAMem")) {
            humanancientpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanancientpool.addModule(new BWAIndex(communicator));
            }
            addBWAMemMapping(humanancientpool);
        }

        if (communicator.isRun_mapping() && communicator.getMapper_to_use().equals("CircularMapper")) {
            humanancientpool.addPredecessor(preprocesspool);
            if (prq.checkForIndices("whole")) {
                humanancientpool.addModule(new BWAIndex(communicator));
            }
            addCircularMapping(humanancientpool);
        }

        if (communicator.isRun_complexityestimation()) {
            addComplexityEstimation(humanancientpool);
        }

        if (communicator.isRmdup_run() && !communicator.isMarkdup_run()) {
            humanancientpool.addModule(new DeDup(communicator));
        }

        if (communicator.isMarkdup_run()) {
            humanancientpool.addModule(new MarkDuplicates(communicator));
        }

        if(communicator.isRun_mapping() && (communicator.isMarkdup_run() || communicator.isRmdup_run())){
            humanancientpool.addModule(new SamtoolsIndex(communicator, SamtoolsIndex.DEDUP));
        }


        if(communicator.isSchmutzi_run()){
            addContaminationEstimation(humanancientpool);
        }


        if (communicator.isRun_coveragecalc()) {
            if(communicator.isSnpcapturedata() || communicator.isRun_mt_capture_mode()){
                humanancientpool.addModule(new QualiMap(communicator,QualiMap.CAPTURE));
            } else {
                humanancientpool.addModule(new QualiMap(communicator));
            }
            humanancientpool.addModule(new MTToNucRatioCalculator(communicator));
        }

        if (communicator.isRun_mapdamage()) {
            humanancientpool.addModule(new MapDamage(communicator));
        }

        if (communicator.isRun_gatksnpcalling() && !communicator.getGatk_caller().equals("ANGSD")) {
            gatkpool = createGATKSNPCallingPipeline();
            gatkpool.addPredecessor(humanancientpool);
            gatkpool.setCurrentFilePath(humanancientpool.getCurrentFilePath());
        }

        if(communicator.isRun_gatksnpcalling() && communicator.getGatk_caller().equals("ANGSD")){
            gatkpool = createANGSDCallingPipeline();
            gatkpool.addPredecessor(humanancientpool);
            gatkpool.setCurrentFilePath(humanancientpool.getCurrentFilePath());
        }

        if (communicator.isRun_gatksnpfiltering()) {
            gatkpool.addModule(new GATKVariantFilter(communicator));
            if (!communicator.isRun_vcf2draft()) {
                gatkpool.addModule(new BGZip(communicator));
                gatkpool.addModule(new Tabix(communicator));
            }
        }

        if (communicator.isRun_vcf2draft()) {
            gatkpool.addModule(new VCF2Genome(communicator));
        }

        if (communicator.isRun_cleanup()) {
            gatkpool.addModule(new CleanUpRedundantData(communicator));
        }

        if (communicator.isRun_reportgenerator()) {
            reportpool.addModule(new ReportGenerator(communicator));
        }

        pools.add(preprocesspool);
        if ( humanancientpool.getModules().size() > 0 ) {
            humanancientpool.addPredecessor(preprocesspool);
            pools.add(humanancientpool);
            reportpool.addPredecessor(humanancientpool);
        } else {
            System.out.println ("Skipping HumanAncient Modules");
        }
        if ( gatkpool.getModules().size() > 0 ) {
            gatkpool.addPredecessor(humanancientpool);
            pools.add(gatkpool);
            reportpool.addPredecessor(gatkpool);
        } else {
            System.out.println ("Skipping GATK Modules");
        }
        if (reportpool.getModules().size() > 0) {
            reportpool.addModule(new ReportGenerator(communicator));
            pools.add(reportpool);
        }
    }

    /**
     * Abstracted methods that are used in all types of pipelines and therefore have been stored in
     * methods.
     */

    /**
     * Adds the BWA mapping modules to the pipeline. First, we determine if we have the following cases:
     * -- Paired End or single ended data
     * -- to be merged or not merged data
     * --> depending on this we decide if sampe or samse have to be used.
     */
    private void addBWAMapping(ModulePool pooltoadd) {
        if (communicator.getMerge_type().equals("PAIRED") && communicator.isRun_clipandmerge() && communicator.isMerge_only_clipping()) {
            pooltoadd.addModule(new BWAAlign(communicator, BWAAlign.PAIREDENDWITHOUTMERGING));
            pooltoadd.addModule(new BWASampe(communicator));
        } else {
            pooltoadd.addModule(new BWAAlign(communicator));
            pooltoadd.addModule(new BWASamse(communicator));
        }

        if (!this.communicator.getMapper_mapquality_filter().equals("0")) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.FILTERED));
            pooltoadd.addModule(new Flagstat(communicator,Flagstat.FILTERED));
        } else {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));
        }
        pooltoadd.addModule(new SamtoolsSort(communicator));


        if (communicator.isRun_mapping_extractmappedandunmapped()) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTMAPPED));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTUNMAPPED));
        }


        pooltoadd.addModule(new SamtoolsIndex(communicator));

    }

    /**
     * Adds the Stampy Mapping modules to the pipeline
     */

    private void addStampyMapping(ModulePool pooltoadd) {
        pooltoadd.addModule(new Stampy(communicator));

        if (!this.communicator.getMapper_mapquality_filter().equals("0")) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.FILTERED));
            pooltoadd.addModule(new Flagstat(communicator));

        } else {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));

        }

        pooltoadd.addModule(new SamtoolsSort(communicator));


        pooltoadd.addModule(new AddOrReplaceReadGroups(communicator)); //Required, Stampy does not add RG information at all -.-
        if (communicator.isRun_mapping_extractmappedandunmapped()) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTMAPPED));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTUNMAPPED));
        }

    }

    /**
     * Adds the CircularMapping Modules
     *
     * @param pooltoadd
     */

    private void addCircularMapping(ModulePool pooltoadd) {
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        if (prq.checkForIndices("CircGenerator")) {
            pooltoadd.addModule(new CircularMapperGenerator(communicator));
        }
        if (prq.checkForIndices("CircularMapper")) {
            pooltoadd.addModule(new BWAIndex(communicator, BWAIndex.MT));
        }
        pooltoadd.addModule(new BWAAlign(communicator, BWAAlign.MT));
        pooltoadd.addModule(new BWASamse(communicator, BWASamse.SAMSEMT));
        if(communicator.isRun_mt_capture_mode()){
            pooltoadd.addModule(new CircularMapperRealigner(communicator,CircularMapperRealigner.FILTERED));
        } else {
            pooltoadd.addModule(new CircularMapperRealigner(communicator, CircularMapperRealigner.DEFAULT));
        }

        if (this.communicator.getMapper_mapquality_filter() != "0") {
            pooltoadd.addModule(new Flagstat(communicator));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.FILTERED));
            pooltoadd.addModule(new Flagstat(communicator));
        } else {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));
        }
        pooltoadd.addModule(new SamtoolsSort(communicator));

        if (communicator.isRun_mapping_extractmappedandunmapped()) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTMAPPED));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTUNMAPPED));
        }
    }

    /**
     * Adds the BWA Mem Mapping modules to the pipeline
     */

    private void addBWAMemMapping(ModulePool pooltoadd) {
        if (communicator.isPairmenttype() && !communicator.isRun_clipandmerge()) {
            pooltoadd.addModule(new BWAMem(communicator, BWAMem.PAIREDENDWITHOUTMERGE));
        } else {
            pooltoadd.addModule(new BWAMem(communicator));
        }

        if (!this.communicator.getMapper_mapquality_filter().equals("0")) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.FILTERED));
            pooltoadd.addModule(new Flagstat(communicator));

        } else {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.DEFAULT));
            pooltoadd.addModule(new Flagstat(communicator));

        }

        pooltoadd.addModule(new SamtoolsSort(communicator));
        pooltoadd.addModule(new AddOrReplaceReadGroups(communicator));

        if (communicator.isRun_mapping_extractmappedandunmapped()) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTMAPPED));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTUNMAPPED));
        }


    }

    /**
     * Adds the BT2 Mapping modules to the pipeline
     */

    private void addBT2Mapping(ModulePool pooltoadd) {
        pooltoadd.addModule(new Bowtie2(communicator));
        pooltoadd.addModule(new SamtoolsView(communicator));
        pooltoadd.addModule(new AddOrReplaceReadGroups(communicator));
        pooltoadd.addModule(new Flagstat(communicator));

        if (communicator.isRun_mapping_extractmappedandunmapped()) {
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTMAPPED));
            pooltoadd.addModule(new SamtoolsView(communicator, SamtoolsView.EXTRACTUNMAPPED));
        }
        pooltoadd.addModule(new SamtoolsSort(communicator));
    }

    /**
     * Adds contamination estimation to the pipeline
     */

    private void addContaminationEstimation(ModulePool pooltoadd){
        pooltoadd.addModule(new SamtoolsFillmd(communicator));
        pooltoadd.addModule(new SamtoolsIndex(communicator,SamtoolsIndex.SCHMUTZI));
        if(communicator.isRun_mt_capture_mode()){
            pooltoadd.addModule(new ContaminationEstimator(communicator));
            pooltoadd.addModule(new ContaminationEstimatorMT(communicator,ContaminationEstimatorMT.DEFAULT));
            pooltoadd.addModule(new ContaminationEstimatorMT(communicator, ContaminationEstimatorMT.NOTPREDC));
        } else {
            pooltoadd.addModule(new ContaminationEstimator(communicator));
        }
    }


    /**
     * Adds the ANGSD low coverage genotype likelihoods generation pipeline
     */

    private ModulePool createANGSDCallingPipeline(){
        ModulePool mp = new ModulePool();

        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        if (prq.checkForIndices("SeqDict")) {
            mp.addModule(new CreateSequenceDictionary(communicator));
        }
        if (prq.checkForIndices("faidx")) {
            mp.addModule(new SamtoolFaidx(communicator));
        }
        mp.addModule(new GATKRealignerTargetCreator(communicator));
        mp.addModule(new GATKIndelRealigner(communicator));

        mp.addModule(new BAMListCreator(communicator));
        if(communicator.isAngsd_create_fasta() && communicator.getAngsd_fasta_callmethod().equals("1")){
            mp.addModule(new ANGSDGenotyper(communicator, ANGSDGenotyper.WITHFASTA));
        } else if(communicator.isAngsd_create_fasta()) {
            mp.addModule(new ANGSDGenotyper(communicator, ANGSDGenotyper.WITHFASTACOUNTS));
        } else {
            mp.addModule(new ANGSDGenotyper(communicator));
        }
        return mp;
    }




    /**
     * Adds the GATK SNP Calling Modules to the pipeline
     */

    private ModulePool createGATKSNPCallingPipeline() {

        ModulePool mp = new ModulePool();
        PrerequisitesChecker prq = new PrerequisitesChecker(communicator);
        if (prq.checkForIndices("SeqDict")) {
            mp.addModule(new CreateSequenceDictionary(communicator));
        }
        if (prq.checkForIndices("faidx")) {
            mp.addModule(new SamtoolFaidx(communicator));
        }
        mp.addModule(new GATKRealignerTargetCreator(communicator));
        mp.addModule(new GATKIndelRealigner(communicator));

        if (communicator.getGatk_caller().equals("HaplotypeCaller")) {
            mp.addModule(new GATKHaplotypeCaller(communicator));
        } else {
            if (communicator.isGatk_emit_all_confident_sites()) {
                if (communicator.isDbsnpreference()) {
                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_CONF_SITES_DBSNP));
                } else {
                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_CONF_SITES_NODBSNP));
                }
            } else if (communicator.isGatk_emit_all_sites()) {
                if (communicator.isDbsnpreference()) {
                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_ALL_SITES_DBSNP));

                } else {

                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_ALL_SITES_NODBSNP));
                }
            } else if (!communicator.isGatk_emit_all_confident_sites() && !communicator.isGatk_emit_all_sites())


                if (communicator.isDbsnpreference()) {

                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_VARIANT_DBSNP));
                } else {

                    mp.addModule(new GATKUnifiedGenotyper(communicator, GATKUnifiedGenotyper.EMIT_VARIANT_NODBSNP));
                }
            }


        if (!communicator.isRun_vcf2draft() && !communicator.isRun_gatksnpfiltering() && !communicator.isSnpcapturedata()){
            //then we need to keep the VCF files untouched!
            mp.addModule(new BGZip(communicator));
            mp.addModule(new Tabix(communicator));
        }


    return mp;

}


    private void addComplexityEstimation(ModulePool pooltoadd) {
        pooltoadd.addModule(new PreseqCCurveCalculation(communicator));
        pooltoadd.addModule(new PreseqLCExtrapCalculation(communicator));
    }


}
