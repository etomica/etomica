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
import etomica.data.types.DataGroup;
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
import etomica.species.SpeciesSpheresMono;
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
    public final SpeciesSpheresMono species;
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

        potential = new P2HardSphere(space);
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
        final String APP_NAME = "HSMD";

        HSMDParam params = new HSMDParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doGraphics = true;
        }
        final HSMDWidom sim = new HSMDWidom(params);

        MeterPressureHard meterP = new MeterPressureHard(sim.integrator);

        MeterRDF meterRDF = new MeterRDF(sim.space);
        meterRDF.getXDataSource().setNValues(1000);
        meterRDF.getXDataSource().setXMax(2);
        meterRDF.setBox(sim.box);
        meterRDF.setResetAfterData(true);
        DataFork forkRDF = new DataFork();

        MeterWidomInsertion meterWidom = new MeterWidomInsertion(sim.space, sim.getRandom());
        meterWidom.setIntegrator(sim.integrator);
        meterWidom.setNInsert(params.nAtoms / 5);
        meterWidom.setSpecies(sim.species);

        if (params.doGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1000);

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, 10));

            AccumulatorAverageFixed accRDF = new AccumulatorAverageFixed(10);
            accRDF.setPushInterval(1);
            forkRDF.addDataSink(accRDF);

            DataPumpListener pumpRDF = new DataPumpListener(meterRDF, forkRDF, 1000);
            sim.integrator.getEventManager().addListener(pumpRDF);

            AccumulatorAverageCollapsing accP = new AccumulatorAverageCollapsing(200);
            DataPumpListener pumpP = new DataPumpListener(meterP, accP, 100);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(accP);
            simGraphic.add(displayP);

            AccumulatorAverageCollapsing accWidom = new AccumulatorAverageCollapsing(200);
            DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom, 10);
            DisplayTextBoxesCAE displayWidom = new DisplayTextBoxesCAE();
            displayWidom.setAccumulator(accWidom);
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
            displayMu.setLabel("Chemical Potnetial");
            displayMu.setDoShowCurrent(false);
            simGraphic.add(displayMu);

            simGraphic.makeAndDisplayFrame(APP_NAME);

            simGraphic.getController().getDataStreamPumps().add(pumpRDF);
            simGraphic.getController().getDataStreamPumps().add(pumpP);
            simGraphic.getController().getDataStreamPumps().add(pumpWidom);
        }
    }

    /**
     * Inner class for parameters understood by the HSMDCavity constructor
     */
    public static class HSMDParam extends ParameterBase {
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
