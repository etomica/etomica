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
import etomica.graphics.DisplayPlot;
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

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
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

        MeterRDF meterRDF = new MeterRDF(sim.space);
        meterRDF.getXDataSource().setNValues(1000);
        meterRDF.getXDataSource().setXMax(2);
        meterRDF.setBox(sim.box);
        meterRDF.setResetAfterData(true);
        DataFork forkRDF = new DataFork();

        DataProcessorCavity cavityProcessor = new DataProcessorCavity(sim.integrator);
        MeterCavityMapped meterCavityMapped = new MeterCavityMapped(sim.integrator);
        meterCavityMapped.setResetAfterData(true);

        MeterContactRatio meterCR = new MeterContactRatio(sim.integrator);

        if (params.doGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1000);

            ColorSchemePaired colorScheme = new ColorSchemePaired(sim.potential);
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, 10));
            sim.integrator.addCollisionListener(meterCavityMapped);
            DisplayPlot cavityPlot = new DisplayPlot();
            cavityPlot.getDataSet().setUpdatingOnAnyChange(true);
            cavityPlot.getPlot().setYLog(true);
            forkRDF.addDataSink(cavityProcessor);

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

            DataPumpListener pumpRDF = new DataPumpListener(meterRDF, forkRDF, 100000);
            sim.integrator.getEventManager().addListener(pumpRDF);
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
        }
    }

    /**
     * Inner class for parameters understood by the HSMDCavity constructor
     */
    public static class CavityParam extends ParameterBase {
        /**
         * Number of atoms, default = 256
         */
        public int nAtoms = 256;
        /**
         * Packing fraction, default = 0.35
         */
        public double eta = 0.35;
        /**
         * Flag indicating whether neighbor list is to be used, default = true
         */
        public boolean useNeighborLists = true;

        public boolean doGraphics = false;
    }

}
