/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.math.DoubleRange;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.dimensions.Null;
import etomica.virial.*;
import etomica.virial.overlap.DataProcessorVirialOverlap;
import etomica.virial.overlap.DataVirialOverlap;

import java.io.*;

/**
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 * 
 * @author kofke, Andrew Schultz
 */
public class SimulationVirialOverlap2 extends Simulation {

    /**
     * This constructor will create your simulation class, but you may call
     * set methods before using it.  When you are done calling set methods,
     * you must call init() before using it.
     */
    public SimulationVirialOverlap2(Space aSpace, ISpecies species, int nMolecules,
                                    double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
        this(aSpace, new ISpecies[]{species}, new int[]{nMolecules}, temperature, refCluster, targetCluster);
    }


    /**
     * This constructor will create your simulation class, but you may call
     * set methods before using it.  When you are done calling set methods,
     * you must call init() before using it.
     */
    public SimulationVirialOverlap2(Space aSpace, ISpecies[] species, int[] nMolecules,
                                    double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
        super(aSpace);
        this.species = species;
        this.temperature = temperature;
        this.nMolecules = nMolecules;
        valueClusters = new ClusterAbstract[]{refCluster, targetCluster};
        sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)};
        meters = new MeterVirial[2];
        integrators = new IntegratorMC[2];
        dpVirialOverlap = new DataProcessorVirialOverlap[2];
        box = new BoxCluster[2];
        accumulators = new AccumulatorRatioAverageCovarianceFull[2];
        extraTargetClusters = new ClusterAbstract[0];
        boxFactory = new BoxClusterFactory();
    }


    /*
     * If this constructor is used to instantiate the simulation, then doWiggle is set to false, and 
     * ClusterAbstract[] is set to {refCluster,targetCluster}
     */
    public SimulationVirialOverlap2(Space aSpace, ISpecies species,
                                    double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
        this(aSpace, new ISpecies[]{species}, new int[]{refCluster.pointCount()}, temperature, refCluster, targetCluster);
        init();
    }

    // this constructor allows you to specify doWiggle=true
    public SimulationVirialOverlap2(Space aSpace, ISpecies species,
                                    double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster, boolean doWiggle) {
        this(aSpace,new ISpecies[]{species}, new int[]{refCluster.pointCount()},temperature,refCluster,targetCluster);
        setDoWiggle(doWiggle);
        init();
    }
    
    // this constructor allows you to specify your own sampleClusters
    public SimulationVirialOverlap2(Space aSpace, ISpecies species,
                                    double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters, boolean doWiggle) {
        this(aSpace, new ISpecies[]{species}, new int[]{aValueClusters[0].pointCount()}, temperature, aValueClusters[0], aValueClusters[1]);
        setDoWiggle(doWiggle);
        setSampleClusters(aSampleClusters);
        init();
    }

    // this constructor allows you to perform the calculation for a mixture
    public SimulationVirialOverlap2(Space aSpace, ISpecies[] species, int[] nMolecules,
                                    double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters, boolean doWiggle) {
        this(aSpace, species, nMolecules, temperature, aValueClusters[0], aValueClusters[1]);
        setSampleClusters(aSampleClusters);
        setDoWiggle(doWiggle);
        init();
    }

    // this constructor allows you to perform the calculation for a mixture or a flexible molecule (with an alternate/ghost molecule)
    // this constructor also allows you to specify extra target diagrams to be calculated during the simulation
    public SimulationVirialOverlap2(Space aSpace, ISpecies[] species, int[] nMolecules,
                                    double temperature, final ClusterAbstract[] aValueClusters, final ClusterAbstract[] extraTargetClusters, final ClusterWeight[] aSampleClusters, boolean doWiggle) {
        this(aSpace, species, nMolecules, temperature, aValueClusters[0], aValueClusters[1]);
        setSampleClusters(aSampleClusters);
        setDoWiggle(doWiggle);
        setExtraTargetClusters(extraTargetClusters);
        init();
	}
    
    public void setDoWiggle(boolean newDoWiggle) {
        if (initialized) throw new RuntimeException("too late");
        doWiggle = newDoWiggle;
    }
    
    public boolean getDoWiggle() {
        return doWiggle;
    }

    public void setExtraTargetClusters(ClusterAbstract[] extraTargetClusters) {
        if (initialized) throw new RuntimeException("too late");
        this.extraTargetClusters = extraTargetClusters;
    }
    
    public ClusterAbstract[] getExtraTargetClusters() {
        return extraTargetClusters;
    }
    
    public void setSampleClusters(ClusterWeight[] sampleClusters) {
        if (initialized) throw new RuntimeException("too late");
        this.sampleClusters[0] = sampleClusters[0];
        this.sampleClusters[1] = sampleClusters[1];
    }
    
    public ClusterWeight[] getSampleClusters() {
        return sampleClusters;
    }
    
    public void setBoxFactory(BoxClusterFactory newBoxFactory) {
        if (initialized) throw new RuntimeException("too late");
        boxFactory = newBoxFactory;
    }
    
    public BoxClusterFactory getBoxFactory() {
        return boxFactory;
    }
    
    public void init() {
        if (initialized) throw new RuntimeException("you can only call me once");
        // we aren't actually initialized yet, but we will be unless we crash.
        // if we crash, we shouldn't get called again!
        initialized = true;

        numExtraTargetClusters = extraTargetClusters.length;
        PotentialMaster potentialMaster = new PotentialMaster();
        boolean doRotate = false;
        boolean multiAtomic = false;
        for (int i=0; i<species.length; i++) {
            addSpecies(species[i]);
            if (!(species[i] instanceof SpeciesSpheresMono) || species[i] instanceof SpeciesSpheresRotating ) {
                doRotate = true;
            }
            if (!(species[i] instanceof SpeciesSpheresMono || species[i] instanceof SpeciesSpheresRotating)) {
                multiAtomic = true;
            }
        }
        accumulatorPumps = new DataPumpListener[2];
        mcMoveTranslate = new MCMoveBoxStep[2];
        if (doRotate) {
            mcMoveRotate = new MCMoveBoxStep[2];
        }
        if (doWiggle) {
            mcMoveWiggle = new MCMoveBoxStep[2];
        }
        
        blockSize = 1000;
        
        for (int iBox=0; iBox<2; iBox++) {
            // integrator for iBox samples based on iBox cluster
            box[iBox] = boxFactory.makeBox(space, sampleClusters[iBox]);
            addBox(box[iBox]);
            for (int i = 0; i < species.length; i++) {
                box[iBox].setNMolecules(species[i], nMolecules[i]);
            }

            integrators[iBox] = new IntegratorMC(this, potentialMaster, box[iBox]);
            integrators[iBox].setTemperature(temperature);
            integrators[iBox].getMoveManager().setEquilibrating(true);

            MCMoveManager moveManager = integrators[iBox].getMoveManager();

            if (!multiAtomic) {
                mcMoveTranslate[iBox] = new MCMoveClusterAtomMulti(random, space);
                moveManager.addMCMove(mcMoveTranslate[iBox]);

                if (doRotate) {
                    mcMoveRotate[iBox] = new MCMoveClusterAtomRotateMulti(random, space);
                    moveManager.addMCMove(mcMoveRotate[iBox]);
                }
            } else {
                mcMoveRotate[iBox] = new MCMoveClusterRotateMoleculeMulti(random, space);
                mcMoveRotate[iBox].setStepSize(Math.PI);
                moveManager.addMCMove(mcMoveRotate[iBox]);
                mcMoveTranslate[iBox] = new MCMoveClusterMoleculeMulti(this, space);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
                if (doWiggle) {
                    // we can use the bending move if none of the molecules has more than 3 atoms
                    boolean doBend = true;
                    for (int i = 0; i < species.length; i++) {
                        if (box[iBox].getNMolecules(species[i]) > 0 && box[iBox].getMoleculeList(species[i]).getMolecule(0).getChildList().getAtomCount() > 3) {
                            doBend = false;
                        }
                    }
                    if (doBend) {
                        mcMoveWiggle[iBox] = new MCMoveClusterAngleBend(potentialMaster, random, 0.5, space);
                    } else {
                        mcMoveWiggle[iBox] = new MCMoveClusterWiggleMulti(this, potentialMaster, valueClusters[0].pointCount(), space);
                    }
                    moveManager.addMCMove(mcMoveWiggle[iBox]);
                }
            }

            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.initializeCoordinates(box[iBox]);
            if (iBox == 0) {
                meters[iBox] = new MeterVirial(new ClusterAbstract[]{valueClusters[0], sampleClusters[1].makeCopy()});
            } else {
                ClusterAbstract[] allClustersForTarget = new ClusterAbstract[extraTargetClusters.length + 2];
                allClustersForTarget[0] = valueClusters[1];
                System.arraycopy(extraTargetClusters, 0, allClustersForTarget, 1, extraTargetClusters.length);
                allClustersForTarget[allClustersForTarget.length - 1] = sampleClusters[0].makeCopy();
                meters[iBox] = new MeterVirial(allClustersForTarget);
            }
            meters[iBox].setBox(box[iBox]);
            dpVirialOverlap[iBox] = new DataProcessorVirialOverlap(11, iBox == 0);
            accumulators[iBox] = new AccumulatorRatioAverageCovarianceFull(blockSize);
            dpVirialOverlap[iBox].setDataSink(accumulators[iBox]);
            accumulatorPumps[iBox] = new DataPumpListener(meters[iBox], dpVirialOverlap[iBox]);
            integrators[iBox].getEventManager().addListener(accumulatorPumps[iBox]);
        }
        
        setRefPref(1,5);
        integratorOS = new IntegratorOverlap(integrators);
        integratorOS.setNumSubSteps(1000);
        integratorOS.setEventInterval(1);
        integratorOS.setAggressiveAdjustStepFraction(true);
        ai = new ActivityIntegrate(integratorOS);
        getController().addAction(ai);
        
        dvo = new DataVirialOverlap(dpVirialOverlap[0], accumulators[0], accumulators[1]);
        integratorOS.setReferenceFracSource(dvo);

    }

    public void setAccumulatorBlockSize(long newBlockSize) {
        blockSize = newBlockSize;
        for (int i=0; i<2; i++) {
            accumulators[i].setBlockSize(newBlockSize);
        }
        // reset the integrator so that it will re-adjust step frequency
        // and ensure it will take enough data for both ref and target
        integratorOS.reset();
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        dpVirialOverlap[0].setBennetParam(refPrefCenter,span);
        dpVirialOverlap[1].setBennetParam(refPrefCenter,span);
    }

    /**
     * Sets the number of alpha values used for the production stage of the
     * simulation.  The default value (1) is sufficient most of the time.
     */
    public void setNumAlpha(int newNumAlpha) {
        numAlpha = newNumAlpha;
        dpVirialOverlap[0].setNumAlpha(newNumAlpha);
        dpVirialOverlap[1].setNumAlpha(newNumAlpha);
    }

    /**
     * Returns the number of alpha values used for the production stage of the
     * simulation
     */
    public int getNumAlpha() {
        return numAlpha;
    }

    protected void initBlockAccumulator() {
        blockAccumulator = new AccumulatorAverageCovariance(1, true);
        DataProcessor dpRatio = new DataProcessor() {
            
            DataDoubleArray data = new DataDoubleArray(extraTargetClusters.length+1);
            
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{extraTargetClusters.length+1});
                return dataInfo;
            }
            
            protected IData processData(IData inputData) {
                double[] x = data.getData();
                double oavg = inputData.getValue(x.length);
                for (int i=0; i<x.length; i++) {
                    x[i] = inputData.getValue(i)/oavg;
                }
                return data;
            }
        };
        accumulators[1].setBlockDataSink(dpRatio);
        dpRatio.setDataSink(blockAccumulator);
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        dpVirialOverlap[0].setNumAlpha(1);
        dpVirialOverlap[0].setBennetParam(newRefPref, 1);
        dpVirialOverlap[1].setNumAlpha(1);
        dpVirialOverlap[1].setBennetParam(newRefPref, 1);
        if (extraTargetClusters.length > 0) {
            initBlockAccumulator();
        }
    }

    /**
     * Causes a progress report (coefficient value and uncertainty) 10 times
     * during the course of the simulation.
     * 
     * The listener is returned.
     */
    public IntegratorListener addProgressListener(final double HSB) {
        return addProgressListener(HSB, false);
    }

    /**
     * Causes a progress report (coefficient value and uncertainty) 10 times
     * during the course of the simulation.  With full=true, the full output
     * (as would be printed at the end of the simulation) is printed.
     * 
     * The listener is returned.
     */
    public IntegratorListener addProgressListener(final double HSB, final boolean full) {
        IntegratorListener progressReport = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                long interval = ai.getMaxSteps()/10;
                if (integratorOS.getStepCount() % interval != 0) return;
                System.out.print(integratorOS.getStepCount()+" steps: ");
                if (full) {
                    System.out.println();
                    printResults(HSB);
                }
                else {
                    double[] ratioAndError = dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB+", error: "+ratioAndError[1]*HSB);
                }
            }
            
            public void integratorInitialized(IntegratorEvent e) {}
        };
        integratorOS.getEventManager().addListener(progressReport);
        return progressReport;
    }

    /**
     * Sets up a histogram for the target system that measures how often
     * maximum separation distances are visited.  The histogram is printed
     * out 10 times during the simulation and may also be printed at the
     * end of the simulation via printTargetHistogram().
     */
    public void setupTargetHistogram() {
        targHist = new HistogramSimple(90, new DoubleRange(-1, 8));
        targPiHist = new HistogramNotSoSimple(90, new DoubleRange(-1, 8));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = box[1].getCPairSet();
                int nPoints = box[1].getMoleculeList().getMoleculeCount();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                }
                else {
                    r -= 1;
                }
                double pi = box[1].getSampleCluster().value(box[1]);
                targHist.addValue(r);
                targPiHist.addValue(r, pi);
            }

            public void integratorInitialized(IntegratorEvent e) {}
        };

        IntegratorListener histReport = new IntegratorListener() {
            public void integratorInitialized(IntegratorEvent e) {}
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorStepFinished(IntegratorEvent e) {
                long interval = ai.getMaxSteps()/10;
                if (integratorOS.getStepCount() % interval != 0) return;
                printTargetHistogram();
            }
        };
        integrators[1].getEventManager().addListener(histListenerTarget);
        integratorOS.getEventManager().addListener(histReport);
    }

    /**
     * Prints out the target histogram of maximum separation distances as
     * set up by setupTargetHistogram().
     */
    public void printTargetHistogram() {
        if (targHist == null) {
            throw new RuntimeException("you can't print a histogram you don't have");
        }
        System.out.println("**** target histograms ****");
        double[] xValues = targHist.xValues();
        double[] h = targHist.getHistogram();
        double[] hPi = targPiHist.getHistogram();
        for (int i=0; i<xValues.length; i++) {
            if (h[i] > 0) {
                double r = xValues[i];
                double y = h[i];
                double pi = hPi[i];
                if (r < 0) r += 1;
                else {
                    r = Math.exp(r);
                    y /= r;
                }
                System.out.println(r+" "+y+" "+pi);
            }
        }
    }

    public void initRefPref(String fileName, long initSteps) {
        // use the old refpref value as a starting point so that an initial
        // guess can be provided
        double oldRefPref = refPref;
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try { 
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+refPref);
                dpVirialOverlap[0].setNumAlpha(numAlpha);
                dpVirialOverlap[1].setNumAlpha(numAlpha);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            for (int i=0; i<2; i++) {
                integrators[i].getMoveManager().setEquilibrating(true);
            }

            long oldBlockSize = blockSize;
            // 1000 blocks
            long newBlockSize = initSteps*integratorOS.getNumSubSteps()/1000;
            if (newBlockSize < 1000) {
                // make block size at least 1000, even if it means fewer blocks
                newBlockSize = 1000;
            }
            if (newBlockSize > 1000000) {
                // needs to be an int.  1e6 steps/block is a bit crazy.
                newBlockSize = 1000000;
            }
            setAccumulatorBlockSize(newBlockSize);
            dpVirialOverlap[0].setNumAlpha(21);
            dpVirialOverlap[1].setNumAlpha(21);
            setRefPref(oldRefPref,30);
            boolean adjustable = integratorOS.isAdjustStepFraction();
            if (adjustable) {
                // we do this initialization to
                // 1. find alpha
                // 2. get molecules out of their starting configuration
                // 3. find optimal mc move step sizes
                // all of these are about as hard in the reference as in the target system
                // so force integratorOS to run both systems equally.
                integratorOS.setRefStepFraction(0.5);
                integratorOS.setAdjustStepFraction(false);
            }
            ai.setMaxSteps(initSteps);
            ai.actionPerformed();
            if (adjustable) {
                integratorOS.setAdjustStepFraction(true);
            }

            refPref = dvo.getOverlapAverage();
            System.out.println("setting initial ref pref to "+refPref);
            if (Double.isInfinite(refPref) || Double.isNaN(refPref)) {
                dvo.getOverlapAverage();
                throw new RuntimeException("oops");
            }
            setAccumulatorBlockSize(oldBlockSize);
            dpVirialOverlap[0].setNumAlpha(15);
            dpVirialOverlap[1].setNumAlpha(15);
            setRefPref(refPref,4);
            for (int i=0; i<2; i++) {
                integrators[i].reset();
            }
            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        ai.setMaxSteps(initSteps);
        long oldBlockSize = blockSize;
        // 1000 blocks
        long newBlockSize = initSteps*integratorOS.getNumSubSteps()/1000;
        if (newBlockSize < 1000) {
            // make block size at least 1000, even if it means fewer blocks
            newBlockSize = 1000;
        }
        if (newBlockSize > 1000000) {
            // needs to be an int.  1e6 steps/block is a bit crazy.
            newBlockSize = 1000000;
        }
        setAccumulatorBlockSize((int)newBlockSize);
        for (int i=0; i<2; i++) {
            integrators[i].getMoveManager().setEquilibrating(true);
        }
        boolean adjustable = integratorOS.isAdjustStepFraction();
        if (adjustable) {
            // we do this initialization to
            // 1. find alpha
            // 2. get molecules out of their starting configuration
            // 3. find optimal mc move step sizes
            // all of these are about as hard in the reference as in the target system
            // so force integratorOS to run both systems equally.
            integratorOS.setRefStepFraction(0.5);
            integratorOS.setAdjustStepFraction(false);
        }
        ai.actionPerformed();
        if (adjustable) {
            integratorOS.setAdjustStepFraction(true);
        }

        if (refPref == -1) {
            refPref = dvo.getOverlapAverage();
            System.out.println("setting ref pref to "+refPref);
            if (Double.isInfinite(refPref) || Double.isNaN(refPref)) {
                dvo.getOverlapAverage();
                throw new RuntimeException("oops");
            }
            dpVirialOverlap[0].setNumAlpha(numAlpha);
            dpVirialOverlap[1].setNumAlpha(numAlpha);
            setRefPref(refPref,1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dvo.reset();
        }
        setAccumulatorBlockSize(oldBlockSize);
        for (int i=0; i<2; i++) {
            integrators[i].getMoveManager().setEquilibrating(false);
        }
        if (extraTargetClusters.length > 0) {
            initBlockAccumulator();
        }
    }

    public void printResults(double refIntegral) {
        printResults(refIntegral, null);
    }

    public void printResults(double refIntegral, String[] extraNames) {
        double[] ratioAndError = dvo.getAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+" error: "+error);
        System.out.println("abs average: "+ratio*refIntegral+" error: "+error*Math.abs(refIntegral));

        double[] alphaData = dvo.getOverlapAverageAndErrorForAlpha(dvo.getAlphaSource().getAlpha(0));
        System.out.println(String.format("overlap ratio: % 20.15e error: %10.5e", alphaData[0], alphaData[1]));

        DataGroup allYourBase = (DataGroup)accumulators[0].getData();
        IData ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
        IData ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
        IData averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        IData stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        IData correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        IData covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(3));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("reference ratio average: % 20.15e error:  %10.5e  cor: %17.15f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
        System.out.print(String.format("reference average: % 20.15e stdev: %9.4e error: %10.5e cor: %17.15f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("reference overlap average: % 20.15e stdev: %9.4e error: %10.5e cor: % 17.15f\n",
                              averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));
        double refRatioAvg = ratioData.getValue(1);
        double refRatioErr = ratioErrorData.getValue(1);
        
        allYourBase = (DataGroup)accumulators[1].getData();
        ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
        ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
        averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
        stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
        errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
        correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        int n = numExtraTargetClusters;
        correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("target ratio average: % 20.15e  error: %10.5e  cor: % 17.15f\n", ratioData.getValue(n + 1), ratioErrorData.getValue(n + 1), correlationCoef));
        System.out.print(String.format("target average: % 20.15e stdev: %9.4e error: %10.5e cor: % 17.15f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("target overlap average: % 20.15e stdev: %9.4e error: %10.5e cor: % 17.15f\n",
                              averageData.getValue(n+1), stdevData.getValue(n+1), errorData.getValue(n+1), correlationData.getValue(n+1)));

        int nTotal = n+2;
        double oVar = covarianceData.getValue(nTotal*nTotal-1);
        double ed = refRatioErr/refRatioAvg;
        double[] e = new double[n+2];
        double[] ro = new double[n+2];
        double[] var = new double[n+2];
        double[] ocor = new double[n+2];
        double[] dcor = new double[n+1];
        double[] corcoef = new double[n+1];
        double[] rd = new double[n+2];
        for (int i=0; i<n+2; i++) {
            e[i] = errorData.getValue(i)/averageData.getValue(i);
            ro[i] = ratioData.getValue((i*nTotal)+n+1)/ratioErrorData.getValue((i*(n+2))+n+1);
            var[i]= covarianceData.getValue((i)*nTotal+(i));
            ocor[i] = var[i]*oVar == 0 ? 0 : covarianceData.getValue(nTotal*(i)+nTotal-1)/Math.sqrt(var[i]*oVar);
            rd[i] = 1/Math.sqrt(ed*ed + 1/(ro[i]*ro[i]));
        }
        for (int i=1; i<n+1; i++) {
            String name = extraNames == null ? ("Extra " + (i)) : extraNames[i - 1];
            // average is vi/|v| average, error is the uncertainty on that average
            // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)           
            System.out.print(String.format("%s average: % 20.15e  error: %10.15e  ocor: % 17.15f", name, averageData.getValue(i), errorData.getValue(i), ocor[i]));
            System.out.print("  dcor:");
            for (int j=0; j<n+1; j++) {
//                if (i==j) continue;                
                dcor[j] = var[i]*var[j] == 0 ? 0 : covarianceData.getValue((i)*nTotal+(j))/Math.sqrt(var[i]*var[j]);
                System.out.print(String.format(" %20.18f", dcor[j]));
            }
            System.out.println();
            int k = (i)*nTotal+(n+1);
            System.out.print(String.format("%s ratio average: % 20.15e  error: %10.15e  tcor:", name, ratioData.getValue(k), ratioErrorData.getValue(k)));
            
            for (int j=0; j<n+1; j++) {                
                corcoef[j] = ro[i]*ro[j]*(e[n+1]*e[n+1] + e[i]*e[j]*dcor[j] - e[i]*e[n+1]*ocor[i] - e[j]*e[n+1]*ocor[j]);

                System.out.print(String.format(" % 17.15f", corcoef[j]));
            }
            System.out.println();
            
            double avg = ratioData.getValue(k)/refRatioAvg;
            double err = Math.abs(avg)/rd[i];
            System.out.print(String.format("%s full average: %20.15e  error: %10.15e  tcor:", name, refIntegral*avg, Math.abs(refIntegral)*err));

            for ( int j=0; j<n+1;j++){
                int kk = (j) * nTotal + (n + 1);
                double avgj = ratioData.getValue(kk) / refRatioAvg;
                double corrcoeff = Math.signum(avg) * Math.signum(avgj) * rd[i] * rd[j] * (ed * ed + corcoef[j] / (ro[i] * ro[j]));
                
                System.out.print(String.format(" %20.18f", corrcoeff));
            }
            System.out.println();
        }
    }

    protected final ISpecies[] species;
    protected final double temperature;
    protected final ClusterAbstract[] valueClusters;
    protected boolean initialized;
    protected boolean doWiggle;
    protected ClusterAbstract[] extraTargetClusters;
    protected final int[] nMolecules;
    protected BoxClusterFactory boxFactory;
    
	public DataVirialOverlap dvo;
    public final AccumulatorRatioAverageCovarianceFull[] accumulators;
    public AccumulatorAverageCovariance blockAccumulator;
    protected DataPumpListener[] accumulatorPumps;
	protected final ClusterWeight[] sampleClusters;
    public final BoxCluster[] box;
    public final IntegratorMC[] integrators;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveTranslate;
    public MCMoveBoxStep[] mcMoveWiggle;
    public final MeterVirial[] meters;
    public int numExtraTargetClusters;
    public final DataProcessorVirialOverlap[] dpVirialOverlap;
    public ActivityIntegrate ai;
    public IntegratorOverlap integratorOS;
    public double refPref;
    protected long blockSize;
    protected int numAlpha = 1;
    protected HistogramSimple targHist;
    protected HistogramNotSoSimple targPiHist;
    
    public static class BoxClusterFactory {
        public BoxCluster makeBox(Space space, ClusterWeight sampleCluster) {
            return new BoxCluster(sampleCluster, space);
        }
    }
}
