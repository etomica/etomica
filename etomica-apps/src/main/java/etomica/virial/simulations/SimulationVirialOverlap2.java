/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.math.DoubleRange;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.NeighborManagerIntra;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.dimensions.Null;
import etomica.virial.BoxCluster;
import etomica.virial.ConfigurationCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MeterVirial;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterWeight;
import etomica.virial.cluster.ClusterWeightAbs;
import etomica.virial.mcmove.*;
import etomica.virial.overlap.DataProcessorVirialOverlap;
import etomica.virial.overlap.DataVirialOverlap;

import java.io.*;
import java.util.Arrays;

/**
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 * 
 * @author kofke, Andrew Schultz
 */
public class SimulationVirialOverlap2 extends Simulation {

    public final AccumulatorRatioAverageCovarianceFull[] accumulators;
    public final BoxCluster[] box;
    public final IntegratorMC[] integrators;
    public final MeterVirial[] meters;
    public final DataProcessorVirialOverlap[] dpVirialOverlap;
    protected final ISpecies[] species;
    protected final double temperature;
    protected final ClusterAbstract[] valueClusters;
    protected final int[] nMolecules;
    protected final ClusterWeight[] sampleClusters;
    public DataVirialOverlap dvo;
    public AccumulatorAverageCovariance blockAccumulator;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveTranslate;
    public MCMoveBoxStep[] mcMoveWiggle;
    public int numExtraTargetClusters;

    public IntegratorOverlap integratorOS;
    public double refPref;
    protected boolean initialized;
    protected boolean doWiggle;
    protected ClusterAbstract[] extraTargetClusters;
    public DataPumpListener[] accumulatorPumps;
    protected long blockSize;
    protected int numAlpha = 1;
    protected HistogramSimple targHist;
    protected HistogramNotSoSimple targPiHist;
    protected double[] boxLengths = new double[]{0, 0};
    protected PotentialMasterBonding.FullBondingInfo bondingInfo;
    protected Potential2Soft[][] pairPotentials;

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
    public SimulationVirialOverlap2(Space aSpace, SpeciesManager sm, int nMolecules,
                                    double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
        super(aSpace, sm);
        this.temperature = temperature;
        this.species = new ISpecies[0];
        this.nMolecules = new int[]{nMolecules};
        valueClusters = new ClusterAbstract[]{refCluster, targetCluster};
        sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), ClusterWeightAbs.makeWeightCluster(targetCluster)};
        meters = new MeterVirial[2];
        integrators = new IntegratorMC[2];
        dpVirialOverlap = new DataProcessorVirialOverlap[2];
        box = new BoxCluster[2];
        accumulators = new AccumulatorRatioAverageCovarianceFull[2];
        extraTargetClusters = new ClusterAbstract[0];
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
    }

    public boolean getDoWiggle() {
        return doWiggle;
    }

    public void setDoWiggle(boolean newDoWiggle) {
        if (initialized) throw new RuntimeException("too late");
        doWiggle = newDoWiggle;
    }

    public ClusterAbstract[] getExtraTargetClusters() {
        return extraTargetClusters;
    }

    public void setExtraTargetClusters(ClusterAbstract[] extraTargetClusters) {
        if (initialized) throw new RuntimeException("too late");
        this.extraTargetClusters = extraTargetClusters;
    }

    public ClusterWeight[] getSampleClusters() {
        return sampleClusters;
    }

    public void setSampleClusters(ClusterWeight[] sampleClusters) {
        if (initialized) throw new RuntimeException("too late");
        this.sampleClusters[0] = sampleClusters[0];
        this.sampleClusters[1] = sampleClusters[1];
    }

    public void setBoxLengths(double refLength, double targetLength) {
        if (initialized) throw new RuntimeException("too late");
        boxLengths[0] = refLength;
        boxLengths[1] = targetLength;
    }

    public void setIntraPairPotentials(Potential2Soft[][] pairPotentials) {
        if (initialized) throw new RuntimeException("too late");
        this.pairPotentials = pairPotentials;
    }

    public void setBondingInfo(PotentialMasterBonding.FullBondingInfo bondingInfo) {
        this.bondingInfo = bondingInfo;
    }

    public void init() {
        if (initialized) throw new RuntimeException("you can only call me once");
        // we aren't actually initialized yet, but we will be unless we crash.
        // if we crash, we shouldn't get called again!
        initialized = true;

        numExtraTargetClusters = extraTargetClusters.length;
        boolean doRotate = false;
        boolean multiAtomic = false;
        for (int i = 0; i < species.length; i++) {
            addSpecies(species[i]);
        }
        for (ISpecies sp : getSpeciesList()) {
            if (sp.getLeafAtomCount() == 1 && sp.getLeafType() instanceof AtomTypeOriented) {
                doRotate = true;
            }
            if (sp.getLeafAtomCount() > 1) {
                multiAtomic = true;
            }
        }
        accumulatorPumps = new DataPumpListener[2];
        mcMoveTranslate = new MCMoveBoxStep[2];
        if (doRotate || multiAtomic) {
            mcMoveRotate = new MCMoveBoxStep[2];
        }
        if (doWiggle) {
            mcMoveWiggle = new MCMoveBoxStep[2];
        }

        blockSize = 1000;

        for (int iBox=0; iBox<2; iBox++) {
            // integrator for iBox samples based on iBox cluster
            box[iBox] = new BoxCluster(sampleClusters[iBox], space, boxLengths[iBox]);
            addBox(box[iBox]);
            for (ISpecies sp : getSpeciesList()) {
                box[iBox].setNMolecules(sp, nMolecules[sp.getIndex()]);
            }

            PotentialCompute pc = null;
            if (pairPotentials != null) {
                PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box[iBox], bondingInfo);
                PotentialComputePair pcPair = new PotentialComputePair(getSpeciesManager(), box[iBox], new NeighborManagerIntra(box[iBox], bondingInfo), pairPotentials);
                pc = new PotentialComputeAggregate(pmBonding, pcPair);
            }
            else if (bondingInfo != null){
                pc = new PotentialMasterBonding(getSpeciesManager(), box[iBox], bondingInfo);
            }
            else {
                pc = new PotentialComputeAggregate();
            }
            integrators[iBox] = new IntegratorMC(pc, getRandom(), temperature, box[iBox]);
            integrators[iBox].getMoveManager().setEquilibrating(true);

            MCMoveManager moveManager = integrators[iBox].getMoveManager();

            if (!multiAtomic) {
                mcMoveTranslate[iBox] = new MCMoveClusterAtomMulti(random, box[iBox]);
                moveManager.addMCMove(mcMoveTranslate[iBox]);

                if (doRotate) {
                    mcMoveRotate[iBox] = new MCMoveClusterAtomRotateMulti(random, box[iBox]);
                    moveManager.addMCMove(mcMoveRotate[iBox]);
                }
            } else {
                mcMoveRotate[iBox] = new MCMoveClusterRotateMoleculeMulti(random, box[iBox]);
                mcMoveRotate[iBox].setStepSize(Math.PI);
                moveManager.addMCMove(mcMoveRotate[iBox]);
                mcMoveTranslate[iBox] = new MCMoveClusterMoleculeMulti(random, box[iBox]);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
                if (doWiggle) {
                    // we can use the bending move if none of the molecules has more than 3 atoms
                    boolean doBend = true;
                    for (ISpecies sp : getSpeciesList()) {
                        if (box[iBox].getNMolecules(sp) > 0 && box[iBox].getMoleculeList(sp).get(0).getChildList().size() > 3) {
                            doBend = false;
                        }
                    }
                    if (doBend) {
                        mcMoveWiggle[iBox] = new MCMoveClusterAngleBend(pc, random, 0.5, space);
                    } else {
                        mcMoveWiggle[iBox] = new MCMoveClusterWiggleMulti(random, pc, box[iBox]);
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

        dvo = new DataVirialOverlap(dpVirialOverlap[0], accumulators[0], accumulators[1]);
        integratorOS.setReferenceFracSource(dvo);

        this.getController().start();
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
     * Returns the number of alpha values used for the production stage of the
     * simulation
     */
    public int getNumAlpha() {
        return numAlpha;
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
        System.out.println("setting ref pref (explicitly) to " + newRefPref);
        refPref = newRefPref;
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
    public IntegratorListener addProgressListener(final double HSB, long interval) {
        return addProgressListener(HSB, interval, false);
    }

    /**
     * Causes a progress report (coefficient value and uncertainty) 10 times
     * during the course of the simulation.  With full=true, the full output
     * (as would be printed at the end of the simulation) is printed.
     *
     * The listener is returned.
     */
    public IntegratorListener addProgressListener(final double HSB, long interval, final boolean full) {
        IntegratorListener progressReport = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
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
    public void setupTargetHistogram(long interval) {
        targHist = new HistogramSimple(90, new DoubleRange(-1, 8));
        targPiHist = new HistogramNotSoSimple(90, new DoubleRange(-1, 8));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = box[1].getCPairSet();
                int nPoints = box[1].getMoleculeList().size();
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
                if (integratorOS.getStepCount() % interval != 0) return;
                printTargetHistogram();
            }
        };
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

    public double[] findUmbrellaWeights(String uWeightsFileName, long uSteps) {
        int n = extraTargetClusters.length;
        double[] uWeights = new double[n];
        boolean uWeightsFileExists = uWeightsFileName != null && new File(uWeightsFileName).exists();
        if (uWeightsFileExists) {
            String line = null;
            try {
                FileReader fr = new FileReader(uWeightsFileName);
                BufferedReader br = new BufferedReader(fr);
                line = br.readLine();
                fr.close();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            String[] bits = line.split(" ");
            if (bits.length != n) {
                throw new RuntimeException("# of umbrella weights found does not match number of clusters");
            }
            for (int i = 0; i < n; i++) {
                uWeights[i] = Double.parseDouble(bits[i]);
            }
            System.out.println("umbrella weights (from file): " + Arrays.toString(uWeights));
        } else {
            integrators[1].reset();
            // equilibrate
            for (long i = 0; i < uSteps / 2; i++) {
                integrators[1].doStep();
            }
            accumulators[1].reset();
            for (int m = 0; m < n; m++) {
                ((ClusterWeightAbs) extraTargetClusters[m]).setDoAbs(true);
            }

            // collect data to determine umbrella weights
            for (long i = 0; i < uSteps; i++) {
                integrators[1].doStep();
            }

            DataGroup allYourBase = (DataGroup) accumulators[1].getData();
            IData averageData = allYourBase.getData(accumulators[1].AVERAGE.index);
            for (int m = 0; m < n; m++) {
                uWeights[m] = averageData.getValue(1) / averageData.getValue(1 + m);
            }
            System.out.println("umbrella weights: " + Arrays.toString(uWeights));

            if (uWeightsFileName != null) {
                try {
                    FileWriter fw = new FileWriter(uWeightsFileName);
                    fw.write(uWeights[0] + "");
                    for (int i = 1; i < n; i++) {
                        fw.write(" " + uWeights[i]);
                    }
                    fw.write("\n");
                    fw.close();
                } catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }

            integrators[1].reset();
        }

        integrators[1].reset();

        for (int m = 0; m < n; m++) {
            ((ClusterWeightAbs) extraTargetClusters[m]).setDoAbs(false);
        }

        return uWeights;
    }

    public void initRefPref(String fileName, long initSteps) {
        initRefPref(fileName, initSteps, true);
    }

    public void initRefPref(String fileName, long initSteps, boolean runBlocking) {
        Activity activityInitRefPref = new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
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
                        setRefPref(refPref, 1);
                        return;
                    } catch (IOException e) {
                        // file not there, which is ok.
                    }
                }

                double initAlphaSpan = 30;
                while (true) {
                    for (int i = 0; i < 2; i++) {
                        integrators[i].getMoveManager().setEquilibrating(true);
                    }

                    long oldBlockSize = blockSize;
                    // 1000 blocks
                    long newBlockSize = initSteps * integratorOS.getNumSubSteps() / 1000;
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
                    setRefPref(oldRefPref, initAlphaSpan);
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

                    for (int i = 0; i < initSteps; i++) {
                        handle.yield(integratorOS::doStep);
                    }

                    if (adjustable) {
                        integratorOS.setAdjustStepFraction(true);
                    }

                    for (int i = 0; i < 2; i++) {
                        integrators[i].reset();
                    }

                    double newRefPref = dvo.getOverlapAverage();
                    if (Double.isInfinite(newRefPref) || Double.isNaN(newRefPref)) {
                        dvo.getOverlapAverage();
                        throw new RuntimeException("oops refpref "+newRefPref);
                    }
                    if (newRefPref > oldRefPref * Math.exp(initAlphaSpan - 0.01) || newRefPref < oldRefPref * Math.exp(-(initAlphaSpan - 0.001))) {
                        System.out.println("guess for ref pref (" + newRefPref + ") is at the edge of the range considered");
                        oldRefPref = newRefPref;
                        continue;
                    }

                    refPref = newRefPref;
                    System.out.println("setting initial ref pref to " + refPref);
                    setAccumulatorBlockSize(oldBlockSize);
                    dpVirialOverlap[0].setNumAlpha(15);
                    dpVirialOverlap[1].setNumAlpha(15);
                    setRefPref(refPref, 4);
                    // set refPref back to -1 so that later on we know that we've been looking for
                    // the appropriate value
                    refPref = -1;
                    break;
                }
            }
        };

        if (runBlocking) {
            this.getController().runActivityBlocking(activityInitRefPref);
        } else {
            this.getController().addActivity(activityInitRefPref);
        }

    }

    protected void readBoxRestart(BufferedReader br, int iBox) throws IOException {
        String[] bits = br.readLine().split(" ");
        mcMoveTranslate[iBox].setStepSize(Double.parseDouble(bits[0]));
        if (mcMoveRotate != null && mcMoveRotate[iBox] != null) {
            mcMoveRotate[iBox].setStepSize(Double.parseDouble(bits[1]));
        }
        if (mcMoveWiggle != null && mcMoveWiggle[iBox] != null) {
            mcMoveWiggle[iBox].setStepSize(Double.parseDouble(bits[1]));
        }
        IAtomList atoms = box[iBox].getLeafList();
        for (IAtom a : atoms) {
            Vector p = a.getPosition();
            bits = br.readLine().split(" ");
            for (int j = 0; j < p.getD(); j++) {
                p.setX(j, Double.parseDouble(bits[j]));
            }
        }
    }

    public long readRestart(String restartFilename, boolean readRef) {
        System.out.println("reading restart file " + restartFilename);
        try {
            FileReader fr = new FileReader(restartFilename + ".sim");
            BufferedReader br = new BufferedReader(fr);
            if (readRef) {
                readBoxRestart(br, 0);
            }
            readBoxRestart(br, 1);
            integratorOS.readStateFromFile(br);
            br.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        integrators[0].getMoveManager().setEquilibrating(false);
        integrators[1].getMoveManager().setEquilibrating(false);
        accumulators[0].readBlockData(restartFilename + ".ref");
        accumulators[1].readBlockData(restartFilename + ".target");
        long readSteps = accumulators[0].getBlockCount() * accumulators[0].getBlockSize()
                + accumulators[1].getBlockCount() * accumulators[1].getBlockSize();
        System.out.println("read data from " + readSteps + " steps");
        new File(restartFilename + ".sim").delete();
        new File(restartFilename + ".ref").delete();
        new File(restartFilename + ".target").delete();
        return readSteps;
    }

    public void equilibrate(String fileName, long initSteps) {
        this.equilibrate(fileName, initSteps, true);
    }

    public void equilibrate(String fileName, long initSteps, boolean runBlocking) {
        Activity activityEquilibrate = new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                // run a short simulation to get reasonable MC Move step sizes and
                // (if needed) narrow in on a reference preference
                long oldBlockSize = blockSize;
                // 1000 blocks
                long newBlockSize = initSteps * integratorOS.getNumSubSteps() / 1000;
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
                for (int i = 0; i < initSteps; i++) {
                    handle.yield(integratorOS::doStep);
                }
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
                for (int i = 0; i < 2; i++) {
                    integrators[i].getMoveManager().setEquilibrating(false);
                }
                if (extraTargetClusters.length > 0) {
                    initBlockAccumulator();
                }
            }
        };

        if (runBlocking) {
            this.getController().runActivityBlocking(activityEquilibrate);
        } else {
            this.getController().addActivity(activityEquilibrate);
        }


    }

    protected void writeBoxRestart(FileWriter fw, int iBox) throws IOException {
        fw.write("" + mcMoveTranslate[iBox].getStepSize());
        if (mcMoveRotate != null && mcMoveRotate[iBox] != null) {
            fw.write(" " + mcMoveRotate[iBox].getStepSize());
        }
        if (mcMoveWiggle != null && mcMoveWiggle[iBox] != null) {
            fw.write(" " + mcMoveWiggle[iBox].getStepSize());
        }
        fw.write("\n");
        IAtomList atoms = box[iBox].getLeafList();
        for (IAtom a : atoms) {
            Vector p = a.getPosition();
            fw.write(p.getX(0) + " " + p.getX(1) + " " + p.getX(2) + "\n");
        }
    }

    public void setWriteRestart(String restartFilename, boolean writeRef) {
        // we write ref blocks even if we're (otherwise) not writing ref
        accumulators[0].setWriteBlocks(restartFilename + ".ref");
        accumulators[1].setWriteBlocks(restartFilename + ".target");

        IntegratorListener restartWriter = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {
            }

            public void integratorStepFinished(IntegratorEvent e) {
                try {
                    FileWriter fw = new FileWriter(restartFilename + ".sim");
                    if (writeRef) {
                        writeBoxRestart(fw, 0);
                    }
                    writeBoxRestart(fw, 1);
                    integratorOS.writeStateToFile(fw);
                    accumulators[0].writeBlockData();
                    accumulators[1].writeBlockData();
                    fw.close();
                } catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        integratorOS.getEventManager().addListener(restartWriter);
    }

    public void printResults(double refIntegral) {
        printResults(refIntegral, null);
    }

    public void printResults(double refIntegral, String[] extraNames) {
        double[] ratioAndError = dvo.getAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: " + ratio + " error: " + error);
        System.out.println("abs average: "+ratio*refIntegral+" error: "+error*Math.abs(refIntegral));

        double[] alphaData = dvo.getOverlapAverageAndErrorForAlpha(dvo.getAlphaSource().getAlpha(0));
        System.out.println(String.format("overlap ratio: % 20.15e error: %10.15e", alphaData[0], alphaData[1]));

        DataGroup allYourBase = (DataGroup)accumulators[0].getData();
        IData ratioData = allYourBase.getData(accumulators[0].RATIO.index);
        IData ratioErrorData = allYourBase.getData(accumulators[0].RATIO_ERROR.index);
        IData averageData = allYourBase.getData(accumulators[0].AVERAGE.index);
        IData stdevData = allYourBase.getData(accumulators[0].STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(accumulators[0].ERROR.index);
        IData correlationData = allYourBase.getData(accumulators[0].BLOCK_CORRELATION.index);
        IData covarianceData = allYourBase.getData(accumulators[0].BLOCK_COVARIANCE.index);
        double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(3));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("reference ratio average: % 20.15e error:  %10.15e  cor: %17.15f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
        System.out.print(String.format("reference average: % 20.15e stdev: %9.4e error: %10.15e cor: %17.15f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("reference overlap average: % 20.15e stdev: %9.4e error: %10.15e cor: % 17.15f\n",
                              averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));
        double refRatioAvg = ratioData.getValue(1);
        double refRatioErr = ratioErrorData.getValue(1);

        allYourBase = (DataGroup)accumulators[1].getData();
        ratioData = allYourBase.getData(accumulators[1].RATIO.index);
        ratioErrorData = allYourBase.getData(accumulators[1].RATIO_ERROR.index);
        averageData = allYourBase.getData(accumulators[1].AVERAGE.index);
        stdevData = allYourBase.getData(accumulators[1].STANDARD_DEVIATION.index);
        errorData = allYourBase.getData(accumulators[1].ERROR.index);
        correlationData = allYourBase.getData(accumulators[1].BLOCK_CORRELATION.index);
        covarianceData = allYourBase.getData(accumulators[1].BLOCK_COVARIANCE.index);
        int n = numExtraTargetClusters;
        correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("target ratio average: % 20.15e  error: %10.15e  cor: % 17.15f\n", ratioData.getValue(n + 1), ratioErrorData.getValue(n + 1), correlationCoef));
        System.out.print(String.format("target average: % 20.15e stdev: %9.4e error: %10.15e cor: % 17.15f\n",
                              averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
        System.out.print(String.format("target overlap average: % 20.15e stdev: %9.4e error: %10.15e cor: % 17.15f\n",
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

            for (int j = 0; j < n + 1; j++) {
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

    @Override
    public Integrator getIntegrator() {
        return this.integratorOS;
    }
}
