/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterRDF;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Hard sphere simulation that allows one pair of atoms to overlap, which
 * allows the cavity function to be measured.
 *
 * @author Andrew Schultz
 */
public class HSMDWidom extends Simulation {

    /**
     * The Box holding the atoms.
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesGeneral species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSphere potential;

    public final PotentialMaster potentialMaster;

    /**
     * Makes a simulation according to the specified parameters.
     *
     * @param params Parameters as defined by the inner class CavityParam
     */
    public HSMDWidom(HSMDParam params) {

        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox();

        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, sigma * neighborRangeFac, space) : new PotentialMasterMonatomic(this);

        int numAtoms = params.nAtoms;

        integrator = new IntegratorHard(random, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.005);

        getController().addActivity(new ActivityIntegrate(integrator));

        potential = new P2HardSphere(space);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.density);
        inflater.actionPerformed();
        if (space.D() == 3) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        } else {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        }

        if (params.useNeighborLists) {
            NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        } else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
    }

    public static void main(String[] args) {
        final String APP_NAME = "HSMD";

        HSMDParam params = new HSMDParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doGraphics = false;
            params.steps = 1000000;
        }
        final HSMDWidom sim = new HSMDWidom(params);

        int nBins = 500;
        double L = sim.box.getBoundary().getBoxSize().getX(0);
        int xMax = (int) (L * 0.5);
        double L3 = L * Math.sqrt(3);
        if (params.mappingCut < L3) L3 = 2 * params.mappingCut;
        int nBinsLong = (int) Math.ceil(nBins * L3 * 0.5);
        double xMaxMap = nBinsLong / ((double) nBins);
        if (params.mappingCut > 0 && params.mappingCut < xMaxMap) xMaxMap = params.mappingCut;

        MeterRDF meterRDF = new MeterRDF(sim.space, true);
        meterRDF.getXDataSource().setNValues(nBins * xMax);
        meterRDF.getXDataSource().setXMax(xMax);
        meterRDF.setBox(sim.box);
        DataFork forkRDF = new DataFork();

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(sim.space, sim.getRandom());
        meterWidom.setIntegrator(sim.integrator);
        meterWidom.setNInsert(params.nAtoms / 5);
        meterWidom.setSpecies(sim.species);

        DataProcessorGContactP dpPContact = new DataProcessorGContactP(sim.box);
        DataProcessorPContactG dpGContact = new DataProcessorPContactG(sim.box);

        MeterWidomCavity meterWC = new MeterWidomCavity(sim.box, sim.getRandom(), sim.potentialMaster);
        meterWC.setSpecies(sim.species);
        meterWC.setNInsert(params.nAtoms / 5);
        double[] r;
        if (params.doGraphics) {
            r = new double[params.nWidomCavity + 1];
            for (int ir = 0; ir < r.length; ir++) r[ir] = ir * (params.widomCavityMax / params.nWidomCavity);
        } else {
            r = new double[params.nWidomCavity];
            for (int ir = 0; ir < r.length; ir++) r[ir] = (ir + 1) * (params.widomCavityMax / params.nWidomCavity);
        }
        meterWC.setInsertionDistances(r);

        if (params.doGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 100);

            DisplayTextBoxesCAE displayPfromGC = null;
            DisplayTextBoxesCAE displayGCMap = null;
            DisplayPlot gPlot = null;
            DataPumpListener pumpRDF = null;
            if (params.doRDF || params.doMappingRDF) {
                gPlot = new DisplayPlot();
                gPlot.setLabel("g(r)");
                simGraphic.add(gPlot);
            }
            DataProcessorExtract0 gCExtractor = null;
            if (params.doRDF) {
                int rdfInterval = (5 * 200 + params.nAtoms - 1) / params.nAtoms;

                AccumulatorAverageFixed accRDF = new AccumulatorAverageFixed(10000 / rdfInterval);
                accRDF.setPushInterval(1);
                forkRDF.addDataSink(accRDF);
                accRDF.addDataSink(gPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accRDF.AVERAGE});

                DataProcessorFit dpFit = new DataProcessorFit("g(r) fit", 100, 3, false, 1, 1.1);
                accRDF.addDataSink(dpFit, new AccumulatorAverage.StatType[]{accRDF.AVERAGE, accRDF.ERROR});
                dpFit.setDataSink(gPlot.getDataSet().makeDataSink());

                gCExtractor = new DataProcessorExtract0("g(sigma)", true);
                gCExtractor.setErrorSource(dpFit);
                dpFit.addDataSink(gCExtractor);
                DisplayTextBoxesCAE gCDisplay = new DisplayTextBoxesCAE();
                gCExtractor.addDataSink(gCDisplay);
                gCDisplay.setDoShowCurrent(false);
                simGraphic.add(gCDisplay);
                if (params.doMappingRDF) {
                    displayGCMap = new DisplayTextBoxesCAE();
                    displayGCMap.setDoShowCurrent(false);
                    simGraphic.add(displayGCMap);
                }
                gCExtractor.addDataSink(dpGContact);
                displayPfromGC = new DisplayTextBoxesCAE();
                dpGContact.addDataSink(displayPfromGC);
                displayPfromGC.setDoShowCurrent(false);
                displayPfromGC.setLabel("P from g(sigma)");

                pumpRDF = new DataPumpListener(meterRDF, forkRDF, rdfInterval);
                sim.integrator.getEventManager().addListener(pumpRDF);
            } else if (params.doMappingRDF) {
                displayGCMap = new DisplayTextBoxesCAE();
                displayGCMap.setDoShowCurrent(false);
                simGraphic.add(displayGCMap);
            }

            MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
            MeterPressureCollisionCount meterPCC = new MeterPressureCollisionCount(sim.integrator);

            AccumulatorAverageCollapsing accP = new AccumulatorAverageCollapsing(200);
            DataPumpListener pumpP = new DataPumpListener(meterP, accP, 100);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(accP);
            displayP.setDoShowCurrent(false);
            displayP.setDoShowCorrelation(true);
            accP.addDataSink(dpPContact, new AccumulatorAverage.StatType[]{accP.MOST_RECENT, accP.AVERAGE, accP.ERROR});
            DisplayTextBoxesCAE displayContact = new DisplayTextBoxesCAE();
            dpPContact.addDataSink(displayContact);
            displayContact.setLabel("g(sigma) from P");
            displayContact.setDoShowCurrent(false);
            simGraphic.add(displayContact);
            simGraphic.add(displayP);
            if (params.doRDF) simGraphic.add(displayPfromGC);
            AccumulatorAverageCollapsing accPCC = new AccumulatorAverageCollapsing(200);
            DataPumpListener pumpPCC = new DataPumpListener(meterPCC, accPCC, 100);
            sim.integrator.getEventManager().addListener(pumpPCC);
            DisplayTextBoxesCAE displayPCC = new DisplayTextBoxesCAE();
            displayPCC.setDoShowCurrent(false);
            displayPCC.setDoShowCorrelation(true);
            displayPCC.setAccumulator(accPCC);
            simGraphic.add(displayPCC);

            DataProcessorExtract0 gMapCExtractor = null;
            if (params.doMappingRDF) {
                MeterRDFMapped meterRDFMapped = new MeterRDFMapped(sim.integrator);
                meterRDFMapped.getXDataSource().setNValues(xMax * nBins + 1);
                meterRDFMapped.getXDataSource().setXMax(xMax);
                meterRDFMapped.reset();
                meterRDFMapped.setResetAfterData(true);

                AccumulatorAverageFixed accRDFMapped = new AccumulatorAverageFixed(1);
                accRDFMapped.setPushInterval(1);
                DataPumpListener pumpRDFMapped = new DataPumpListener(meterRDFMapped, accRDFMapped, 10000);
                accRDFMapped.addDataSink(gPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accRDFMapped.AVERAGE});
                sim.integrator.getEventManager().addListener(pumpRDFMapped);

                DataProcessorErrorBar rdfMappedError = new DataProcessorErrorBar("mapped g(r)+");
                accRDFMapped.addDataSink(rdfMappedError, new AccumulatorAverage.StatType[]{accRDFMapped.AVERAGE, accRDFMapped.ERROR});
                rdfMappedError.setDataSink(gPlot.getDataSet().makeDataSink());

                gMapCExtractor = new DataProcessorExtract0("mapped g(sigma)", true);
                accRDFMapped.addDataSink(gMapCExtractor, new AccumulatorAverage.StatType[]{accRDFMapped.MOST_RECENT, accRDFMapped.AVERAGE, accRDFMapped.ERROR});
                gMapCExtractor.addDataSink(displayGCMap);
                simGraphic.getController().getDataStreamPumps().add(pumpRDFMapped);
            }

            if (params.doRDF) simGraphic.getController().getDataStreamPumps().add(pumpRDF);
            simGraphic.getController().getDataStreamPumps().add(pumpP);

            if (params.doWidom) {
                AccumulatorAverageCollapsing accWidom = new AccumulatorAverageCollapsing(200);
                DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom, 10);
                DisplayTextBoxesCAE displayWidom = new DisplayTextBoxesCAE();
                displayWidom.setAccumulator(accWidom);
                displayWidom.setDoShowCurrent(false);
                displayWidom.setDoShowCorrelation(true);
                simGraphic.add(displayWidom);
                sim.integrator.getEventManager().addListener(pumpWidom);
                DataProcessorForked widomProcessor = new DataProcessorForked() {
                    DataGroup data = new DataGroup(new DataDouble[]{new DataDouble(), new DataDouble(), new DataDouble()});

                    @Override
                    protected IData processData(IData inputData) {
                        ((DataDouble) data.getData(0)).x = Double.NaN;
                        ((DataDouble) data.getData(1)).x = -Math.log(inputData.getValue(0));
                        ((DataDouble) data.getData(2)).x = inputData.getValue(1) / inputData.getValue(0);
                        return data;
                    }

                    @Override
                    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                        DataDouble.DataInfoDouble subDataInfo = new DataDouble.DataInfoDouble("stuff", Null.DIMENSION);
                        dataInfo = new DataGroup.DataInfoGroup("chemical potential", Energy.DIMENSION, new IDataInfo[]{subDataInfo, subDataInfo, subDataInfo});
                        dataInfo.addTag(tag);
                        return dataInfo;
                    }
                };
                accWidom.addDataSink(widomProcessor, new AccumulatorAverage.StatType[]{accWidom.AVERAGE, accWidom.ERROR});
                DisplayTextBoxesCAE displayMu = new DisplayTextBoxesCAE();
                widomProcessor.addDataSink(displayMu);
                displayMu.setLabel("Chemical Potential");
                displayMu.setDoShowCurrent(false);
                simGraphic.add(displayMu);
                simGraphic.getController().getDataStreamPumps().add(pumpWidom);
            }

            if (params.doWidomCavity) {
                AccumulatorAverageFixed accWC = new AccumulatorAverageFixed(200);
                DataPumpListener pumpWC = new DataPumpListener(meterWC, accWC, 10);
                sim.integrator.getEventManager().addListener(pumpWC);
                DisplayPlot yPlot = new DisplayPlot();
                accWC.addDataSink(yPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accWC.AVERAGE});
                DataProcessorErrorBar dpErrorWC = new DataProcessorErrorBar("y(r)+");
                accWC.addDataSink(dpErrorWC, new AccumulatorAverage.StatType[]{accWC.AVERAGE, accWC.ERROR});
                dpErrorWC.setDataSink(yPlot.getDataSet().makeDataSink());
                yPlot.setLabel("y(r)");
                simGraphic.add(yPlot);
                yPlot.getPlot().setYLog(true);

                if (params.doMappingRDF || params.doRDF) {
                    DataProcessorWidomContact dpwc = new DataProcessorWidomContact();
                    accWC.addDataSink(dpwc, new AccumulatorAverage.StatType[]{accWC.AVERAGE, accWC.ERROR});
                    if (params.doMappingRDF) {
                        gMapCExtractor.addDataSink(dpwc.makeGContactSink());
                    } else {
                        gCExtractor.addDataSink(dpwc.makeGContactSink());
                    }
                    DisplayTextBoxesCAE displayWC = new DisplayTextBoxesCAE();
                    displayWC.setDoShowCurrent(false);
                    dpwc.setDataSink(displayWC);
                    simGraphic.add(displayWC);
                }
            }

            simGraphic.makeAndDisplayFrame(APP_NAME);
            return;
        }

        System.out.println("N: " + params.nAtoms);
        System.out.println("steps: " + params.steps);
        System.out.println("density: " + params.density);

        long steps = params.steps;
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        sim.integrator.resetStepCount();

        AccumulatorAverageFixed accRDFMapped = new AccumulatorAverageFixed(1);
        if (params.doMappingRDF) {
            MeterRDFMapped meterRDFMapped = new MeterRDFMapped(sim.integrator);
            meterRDFMapped.foobar = params.doMappingFoobar;
            if (params.mappingCut > 0) meterRDFMapped.setMappingCut(params.mappingCut);
            meterRDFMapped.getXDataSource().setNValues(nBinsLong + 1);
            meterRDFMapped.getXDataSource().setXMax(xMaxMap);
            meterRDFMapped.reset();
            meterRDFMapped.setResetAfterData(true);

            DataPumpListener pumpRDFMapped = new DataPumpListener(meterRDFMapped, null, (int) (steps / 100));
            meterRDFMapped.setRawSink(accRDFMapped);
            sim.integrator.getEventManager().addListener(pumpRDFMapped);
        }

        AccumulatorAverageFixed accRDF = new AccumulatorAverageFixed(1);
        if (params.doRDF) {
            int rdfInterval = (30 * params.nAtoms + 199) / 200;
            if (params.rdfInterval > 0) rdfInterval = params.rdfInterval;
            int stepsPerBlock = (int) (steps / 100);
            while (rdfInterval > stepsPerBlock && (steps % (stepsPerBlock * 2) == 0)) stepsPerBlock *= 2;
            while (stepsPerBlock % rdfInterval != 0) rdfInterval--;
            System.out.println("RDF interval: " + rdfInterval);
            System.out.println("steps per block: " + stepsPerBlock);
            accRDF.setBlockSize(stepsPerBlock / rdfInterval);

            DataPumpListener pumpRDF = new DataPumpListener(meterRDF, accRDF, rdfInterval);
            sim.integrator.getEventManager().addListener(pumpRDF);
        }

        AccumulatorAverageFixed accWidom = null;
        if (params.doWidom) {
            accWidom = new AccumulatorAverageFixed((int) (steps / 1000));
            DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom, 10);
            sim.integrator.getEventManager().addListener(pumpWidom);
        }

        AccumulatorAverageFixed accWC = null;
        if (params.doWidomCavity) {
            accWC = new AccumulatorAverageFixed((int) (steps / 1000));
            DataPumpListener pumpWC = new DataPumpListener(meterWC, accWC, 10);
            sim.integrator.getEventManager().addListener(pumpWC);
        }

        MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
        MeterPressureCollisionCount meterPCC = new MeterPressureCollisionCount(sim.integrator);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(1);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, (int) (steps / 100));
        sim.integrator.getEventManager().addListener(pumpP);

        AccumulatorAverageFixed accPCC = new AccumulatorAverageFixed(1);
        DataPumpListener pumpPCC = new DataPumpListener(meterPCC, accPCC, (int) (steps / 100));
        sim.integrator.getEventManager().addListener(pumpPCC);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.nanoTime();


        double avgP = accP.getData(accP.AVERAGE).getValue(0);
        double errP = accP.getData(accP.ERROR).getValue(0);
        double corP = accP.getData(accP.BLOCK_CORRELATION).getValue(0);
        System.out.println("Pressure: " + avgP + " err: " + errP + " cor: " + corP);

        double avgPCC = accPCC.getData(accPCC.AVERAGE).getValue(0);
        double errPCC = accPCC.getData(accPCC.ERROR).getValue(0);
        double corPCC = accPCC.getData(accPCC.BLOCK_CORRELATION).getValue(0);
        System.out.println("Pressure(CC): " + avgPCC + " err: " + errPCC + " cor: " + corPCC);

        if (params.doRDF) {
            IData rData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) accRDF.getDataInfo()).getSubDataInfo(0)).getXDataSource().getIndependentData(0);
            IData rdfDataAvg = accRDF.getData(accRDF.AVERAGE);
            IData rdfDataErr = accRDF.getData(accRDF.ERROR);
            IData rdfDataCor = accRDF.getData(accRDF.BLOCK_CORRELATION);
            System.out.println("\nRDF");
            for (int i = 0; i < rData.getLength(); i++) {
                System.out.println(rData.getValue(i) + " " + rdfDataAvg.getValue(i) + " " + rdfDataErr.getValue(i) + " " + rdfDataCor.getValue(i));
            }
        }
        if (params.doMappingRDF) {
            IData rData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) accRDFMapped.getDataInfo()).getSubDataInfo(0)).getXDataSource().getIndependentData(0);
            IData rdfDataAvg = accRDFMapped.getData(accRDFMapped.AVERAGE);
            IData rdfDataErr = accRDFMapped.getData(accRDFMapped.ERROR);
            IData rdfDataCor = accRDFMapped.getData(accRDFMapped.BLOCK_CORRELATION);
            System.out.println("\nmapped RDF");
            for (int i = 0; i < rData.getLength(); i++) {
                System.out.println(rData.getValue(i) + " " + rdfDataAvg.getValue(i) + " " + rdfDataErr.getValue(i) + " " + rdfDataCor.getValue(i));
            }
        }
        if (params.doWidom) {
            double avgExp = accWidom.getData(accWidom.AVERAGE).getValue(0);
            double errExp = accWidom.getData(accWidom.ERROR).getValue(0);
            double corExp = accWidom.getData(accWidom.BLOCK_CORRELATION).getValue(0);
            System.out.println("Widom exp: " + avgExp + "  err: " + errExp + "  cor: " + corExp);
            double mu = -Math.log(avgExp);
            double errMu = errExp / avgExp;
            System.out.println("Widom mu: " + mu + "  err: " + errMu);
        }
        if (params.doWidomCavity) {
            IData rData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) accWC.getDataInfo()).getSubDataInfo(0)).getXDataSource().getIndependentData(0);
            IData wcDataAvg = accWC.getData(accWC.AVERAGE);
            IData wcDataErr = accWC.getData(accWC.ERROR);
            IData wcDataCor = accWC.getData(accWC.BLOCK_CORRELATION);
            System.out.println("\nWidom cavity");
            for (int i = 0; i < rData.getLength(); i++) {
                System.out.println(rData.getValue(i) + " " + wcDataAvg.getValue(i) + " " + wcDataErr.getValue(i) + " " + wcDataCor.getValue(i));
            }
        }

        System.out.println("\ntime: " + (t2 - t1) / 1e9);

    }

    /**
     * Inner class for parameters understood by the HSMDCavity constructor
     */
    public static class HSMDParam extends ParameterBase {
        public long steps = 1000000;
        public int nAtoms = 256;
        public double density = 0.6;
        public int mappingCut = 0;
        public boolean useNeighborLists = true;
        public boolean doGraphics = false;
        public boolean doWidom = false;
        public boolean doRDF = false;
        public int rdfInterval = 0;
        public boolean doMappingRDF = false;
        public boolean doMappingFoobar = false;
        public boolean doWidomCavity = false;
        public int nWidomCavity = 10;
        public double widomCavityMax = 1;
    }

}
