/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterRDF;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Hard sphere simulation that allows one pair of atoms to overlap, which
 * allows the cavity function to be measured.
 *
 * @author Andrew Schultz
 */
public class HSMDCavity extends Simulation {

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
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSphereCavity potential;

    public final PotentialMaster potentialMaster;

    public final ActivityIntegrate activityIntegrate;

    /**
     * Makes a simulation according to the specified parameters.
     *
     * @param params Parameters as defined by the inner class CavityParam
     */
    public HSMDCavity(CavityParam params) {

        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        box = this.makeBox();

        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, sigma * neighborRangeFac, space) : new PotentialMasterMonatomic(this);

        int numAtoms = params.nAtoms;

        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.005);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        potential = new P2HardSphereCavity(space);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
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
        final String APP_NAME = "HSMDCavity";

        CavityParam params = new CavityParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doGraphics = true;
        }
        final HSMDCavity sim = new HSMDCavity(params);

        MeterRDF meterCavity = new MeterCavity(sim.space, sim.potential);
        meterCavity.getXDataSource().setNValues(params.nBins);
        meterCavity.getXDataSource().setXMax(1);
        meterCavity.setBox(sim.box);
        meterCavity.setResetAfterData(true);
        DataFork forkCavity = new DataFork();

        MeterRDF meterRDF = new MeterRDF(sim.space);
        meterRDF.getXDataSource().setNValues(2 * params.nBins);
        meterRDF.getXDataSource().setXMax(2);
        meterRDF.setBox(sim.box);
        meterRDF.setResetAfterData(true);
        DataFork forkRDF = new DataFork();

        DataProcessorCavity cavityProcessor = new DataProcessorCavity(sim.integrator);
        MeterCavityMapped meterCavityMapped = new MeterCavityMapped(sim.integrator);
        meterCavityMapped.setResetAfterData(true);
        meterCavityMapped.getXDataSource().setNValues(params.nBins + 1);

        MeterContactRatio meterCR = new MeterContactRatio(sim.integrator);

        if (params.doGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1000);

            ColorSchemePaired colorScheme = new ColorSchemePaired(sim.potential);
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, 10));
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterCavity, 1));
            sim.integrator.addCollisionListener(meterCavityMapped);
            DisplayPlot cavityPlot = new DisplayPlot();
            cavityPlot.getDataSet().setUpdatingOnAnyChange(true);
            cavityPlot.getPlot().setYLog(true);
            forkCavity.addDataSink(cavityProcessor);

            AccumulatorAverageFixed accMapped = new AccumulatorAverageFixed(1);
            accMapped.setPushInterval(1);
            DataPumpListener pumpCavityMapped = new DataPumpListener(meterCavityMapped, accMapped, 100000);
            accMapped.addDataSink(cavityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accMapped.AVERAGE});
            DataProcessor mappedErr = new DataProcessorErrorBar("mapped y(r)+e");
            accMapped.addDataSink(mappedErr, new AccumulatorAverage.StatType[]{accMapped.AVERAGE, accMapped.ERROR});
            mappedErr.setDataSink(cavityPlot.getDataSet().makeDataSink());

            AccumulatorAverageFixed accRDF = new AccumulatorAverageFixed(1);
            accRDF.setPushInterval(1);
            forkRDF.addDataSink(accRDF);
            accRDF.addDataSink(cavityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accRDF.AVERAGE});

            AccumulatorAverageFixed accCavity = new AccumulatorAverageFixed(1);
            accCavity.setPushInterval(1);
            forkCavity.addDataSink(accCavity);
            accCavity.addDataSink(cavityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accCavity.AVERAGE});
            cavityPlot.setLegend(new DataTag[]{accCavity.getTag()}, "g(r<1)");

            AccumulatorAverageFixed accConv = new AccumulatorAverageFixed(1);
            accConv.setPushInterval(1);
            DataFork cavityFork = new DataFork();
            cavityProcessor.setDataSink(cavityFork);
            cavityFork.addDataSink(accConv);
            DataProcessorFit processorFit = new DataProcessorFit("yFit(r)", meterRDF.getXDataSource().getNValues() / 2 - 1, 8, true);
            accConv.addDataSink(processorFit, new AccumulatorAverage.StatType[]{accConv.AVERAGE, accConv.ERROR});
            processorFit.addDataSink(cavityPlot.getDataSet().makeDataSink());
            accConv.addDataSink(cavityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accConv.AVERAGE});
            DataProcessor convErr = new DataProcessorErrorBar("y(r)+e");
            accConv.addDataSink(convErr, new AccumulatorAverage.StatType[]{accConv.AVERAGE, accConv.ERROR});
            convErr.setDataSink(cavityPlot.getDataSet().makeDataSink());

            DataProcessorFit.DataSourceLastOrder dsOrder = processorFit.makeDataSourceLastOrder();
            DisplayTextBox displayOrder = new DisplayTextBox();
            DataPumpListener pumpOrder = new DataPumpListener(dsOrder, displayOrder, 1000);
            sim.integrator.getEventManager().addListener(pumpOrder);
            simGraphic.add(displayOrder);

            DataProcessorFit.DataSourceChi dsChi = processorFit.makeDataSourceChi();
            DisplayTextBox displayChi = new DisplayTextBox();
            DataPumpListener pumpChi = new DataPumpListener(dsChi, displayChi, 1000);
            sim.integrator.getEventManager().addListener(pumpChi);
            simGraphic.add(displayChi);

            DataPumpListener pumpRDF = new DataPumpListener(meterRDF, forkRDF, 100000);
            sim.integrator.getEventManager().addListener(pumpRDF);
            DataPumpListener pumpCavity = new DataPumpListener(meterCavity, forkCavity, 100000);
            sim.integrator.getEventManager().addListener(pumpCavity);
            sim.integrator.getEventManager().addListener(pumpCavityMapped);
            cavityPlot.setLabel("y(r)");

            DisplayPlot y0Plot = new DisplayPlot();
            DataProcessorExtract0 y0Extractor = new DataProcessorExtract0("conv", null);
            processorFit.addDataSink(y0Extractor);
            DataSourceCountTime dsTime = new DataSourceCountTime(sim.integrator);
            AccumulatorHistory y0History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            y0History.setTimeDataSource(dsTime);
            y0Extractor.addDataSink(y0History);
            y0History.addDataSink(y0Plot.getDataSet().makeDataSink());
            y0Plot.setLabel("y(0)");
            DataProcessorExtract0 map0Extractor = new DataProcessorExtract0("mapped", null);
            accMapped.addDataSink(map0Extractor, new AccumulatorAverage.StatType[]{accMapped.AVERAGE});
            AccumulatorHistory map0History = new AccumulatorHistory(new HistoryCollapsingDiscard());
            map0History.setTimeDataSource(dsTime);
            map0Extractor.addDataSink(map0History);
            map0History.addDataSink(y0Plot.getDataSet().makeDataSink());

            AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(1);
            accForce.setPushInterval(1);
            meterCavityMapped.setForceSink(accForce);
            DisplayPlot forcePlot = new DisplayPlot();
            accForce.addDataSink(forcePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accForce.AVERAGE});
            forcePlot.setLabel("force");
            simGraphic.add(forcePlot);
            forcePlot.getPlot().setYLog(true);

            AccumulatorHistory accCR = new AccumulatorHistory(new HistoryCollapsingAverage());
            accCR.setPushInterval(1);
            DataPumpListener pumpCR = new DataPumpListener(meterCR, accCR, 1000);
            sim.integrator.getEventManager().addListener(pumpCR);
            DisplayPlot plotCR = new DisplayPlot();
            accCR.addDataSink(plotCR.getDataSet().makeDataSink());
            plotCR.setLabel("CR");
            simGraphic.add(plotCR);

            simGraphic.add(cavityPlot);
            simGraphic.add(y0Plot);

            simGraphic.makeAndDisplayFrame(APP_NAME);

            simGraphic.getController().getDataStreamPumps().add(pumpRDF);
            simGraphic.getController().getDataStreamPumps().add(pumpCavityMapped);
            simGraphic.getController().getResetAveragesButton().setPostAction(new IAction() {
                @Override
                public void actionPerformed() {
                    accForce.reset();
                }
            });
            return;
        }

        long steps = params.steps;
        sim.activityIntegrate.setMaxSteps(steps / 10);
        sim.getController().actionPerformed();

        if (params.doConv) {
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, 10));
        }
        if (params.doMapping) {
            sim.integrator.addCollisionListener(meterCavityMapped);
        }

        AccumulatorAverageFixed accMapped = new AccumulatorAverageFixed(1);
        DataPumpListener pumpCavityMapped = new DataPumpListener(meterCavityMapped, null, (int) (steps / 100));
        meterCavityMapped.setForceSink(accMapped);

        AccumulatorAverageFixed accRDF = new AccumulatorAverageFixed(1);
        forkRDF.addDataSink(accRDF);

        DataPumpListener pumpRDF = new DataPumpListener(meterRDF, forkRDF, (int) (steps / 100));
        sim.integrator.getEventManager().addListener(pumpRDF);
        sim.integrator.getEventManager().addListener(pumpCavityMapped);

        AccumulatorAverageFixed accCR = new AccumulatorAverageFixed(1);
        DataPumpListener pumpCR = new DataPumpListener(meterCR, accCR, (int) (steps / 100));
        sim.integrator.getEventManager().addListener(pumpCR);

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getController().reset();
        double t1 = System.nanoTime();
        sim.getController().actionPerformed();
        double t2 = System.nanoTime();

        double avgCR = accCR.getData(accCR.AVERAGE).getValue(0);
        double errCR = accCR.getData(accCR.ERROR).getValue(0);
        double corCR = accCR.getData(accCR.BLOCK_CORRELATION).getValue(0);
        System.out.println("contact ratio average: " + avgCR + "  error: " + errCR + " cor: " + corCR);

        if (params.doConv) {
            DataGroup rdfData = (DataGroup) accRDF.getData();
            IData avgRDF = rdfData.getData(accRDF.AVERAGE.index);
            IData errRDF = rdfData.getData(accRDF.ERROR.index);
            IData corRDF = rdfData.getData(accRDF.BLOCK_CORRELATION.index);
            IData rData = ((DataFunction.DataInfoFunction) meterRDF.getDataInfo()).getXDataSource().getIndependentData(0);
            System.out.println("\nRDF");
            for (int i = 0; i < rData.getLength(); i++) {
                System.out.println(rData.getValue(i) + " " + avgRDF.getValue(i) + " " + errRDF.getValue(i) + " " + corRDF.getValue(i));
            }
        }

        if (params.doMapping) {
            DataGroup yMapData = (DataGroup) accMapped.getData();
            IData avgYMapData = yMapData.getData(accMapped.AVERAGE.index);
            IData errYMapData = yMapData.getData(accMapped.ERROR.index);
            IData corYMapData = yMapData.getData(accMapped.BLOCK_CORRELATION.index);
            IData rData = ((DataFunction.DataInfoFunction) meterCavityMapped.getDataInfo()).getXDataSource().getIndependentData(0);
            System.out.println("\nmapped y(r)-y(sigma)");
            for (int i = 0; i < rData.getLength(); i++) {
                System.out.println(rData.getValue(i) + " " + avgYMapData.getValue(i) + " " + errYMapData.getValue(i) + " " + corYMapData.getValue(i));
            }
        }

        System.out.println("time: " + (t2 - t1) / 1e9);
    }

    /**
     * Inner class for parameters understood by the HSMDCavity constructor
     */
    public static class CavityParam extends ParameterBase {
        public int nAtoms = 256;
        public double eta = 0.35;
        public boolean useNeighborLists = true;
        public boolean doGraphics = false;
        public int nBins = 500;
        public long steps = 100000;
        public boolean doConv = true;
        public boolean doMapping = true;
    }

}
