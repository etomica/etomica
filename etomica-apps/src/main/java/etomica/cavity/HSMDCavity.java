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
import etomica.data.meter.MeterRDF;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
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
        integrator.setTimeStep(0.01);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
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
        meterRDF.getXDataSource().setXMax(2);
        meterRDF.setBox(sim.box);
        DataFork forkRDF = new DataFork();
        DataProcessorRDF rdfProcessor = new DataProcessorRDF(sim.potential.getCollisionDiameter());
        forkRDF.addDataSink(rdfProcessor);

        DataProcessorCavity cavityProcessor = new DataProcessorCavity(sim.integrator, sim.potential);
        MeterCavityMapped meterCavityMapped = new MeterCavityMapped(sim.integrator);
        meterCavityMapped.setResetAfterData(true);

        if (params.doGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);

            ColorSchemePaired colorScheme = new ColorSchemePaired(sim.potential);
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, 10));
            sim.integrator.addCollisionListener(meterCavityMapped);
            DisplayPlot cavityPlot = new DisplayPlot();
            cavityPlot.getDataSet().setUpdatingOnAnyChange(true);
            forkRDF.addDataSink(cavityProcessor);
            AccumulatorAverageFixed accMapped = new AccumulatorAverageFixed(1);
            accMapped.setPushInterval(1);
            DataPumpListener pumpCavityMapped = new DataPumpListener(meterCavityMapped, accMapped, 10000);
            accMapped.addDataSink(cavityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{accMapped.AVERAGE});
            DataProcessor mappedErr = new DataProcessor() {
                protected DataFunction data;

                @Override
                protected IData processData(IData inputData) {
                    IData avg = ((DataGroup) inputData).getData(0);
                    IData err = ((DataGroup) inputData).getData(1);
                    double[] y = data.getData();
                    for (int i = 0; i < y.length; i++) {
                        y[i] = avg.getValue(i) + err.getValue(i);
                    }
                    return data;
                }

                @Override
                protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                    inputDataInfo = ((DataGroup.DataInfoGroup) inputDataInfo).getSubDataInfo(0);
                    data = (DataFunction) inputDataInfo.makeData();
                    IDataInfoFactory factory = inputDataInfo.getFactory();
                    factory.setLabel("mapped y(r)+e");
                    dataInfo = factory.makeDataInfo();
                    dataInfo.addTag(tag);
                    return dataInfo;
                }
            };
            accMapped.addDataSink(mappedErr, new AccumulatorAverage.StatType[]{accMapped.AVERAGE, accMapped.ERROR});
            mappedErr.setDataSink(cavityPlot.getDataSet().makeDataSink());

            forkRDF.addDataSink(cavityPlot.getDataSet().makeDataSink());
            cavityProcessor.setDataSink(cavityPlot.getDataSet().makeDataSink());
            DataPumpListener pumpRDF = new DataPumpListener(meterRDF, forkRDF, 1000);
            sim.integrator.getEventManager().addListener(pumpRDF);
            sim.integrator.getEventManager().addListener(pumpCavityMapped);
            cavityPlot.setLabel("cavity");
            simGraphic.add(cavityPlot);
            simGraphic.makeAndDisplayFrame(APP_NAME);

            simGraphic.getController().getResetAveragesButton().setPostAction(new IAction() {
                @Override
                public void actionPerformed() {
                    sim.integrator.resetStepCount();
                    sim.potential.resetInternalCount();
                    meterRDF.reset();
                    meterCavityMapped.reset();
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
