/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.FunctionDifferentiable;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.List;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class MappedSingletDensity extends Simulation {

    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public MappedSingletDensity() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public MappedSingletDensity(SimParams params) {
        super(Space3D.getInstance());

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        double rc = 3;
        potentialMaster = new PotentialMasterCell(this, rc, space);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones(space);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.addPotential(p2, new AtomType[]{atomType, atomType});

        if (params.field == Field.SINE) {
            P1Sine p1 = new P1Sine(space, 5, params.temperature);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        } else if (params.field == Field.PARABOLIC) {
            throw new RuntimeException("not yet");
        }

        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(params.temperature);

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        potentialMaster.setCellRange(2);
        potentialMaster.reset();
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.field = Field.SINE;
            // modify parameters here for interactive testing
        }

        MappedSingletDensity sim = new MappedSingletDensity(params);
        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.activityIntegrate.setMaxSteps(steps / 10);
        sim.activityIntegrate.actionPerformed();
        System.out.println("equilibration finished");

        // data collection
        MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
        densityMeter.setBox(sim.box());
        densityMeter.setProfileDim(2);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        densityMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(densityMeter, acc, interval);
        sim.getIntegrator().getEventManager().addListener(pump);

        MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
        AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, interval);
        sim.getIntegrator().getEventManager().addListener(pumpForce);

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getIntegrator().resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.actionPerformed();

        long t2 = System.currentTimeMillis();

        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();
            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            }
            else {
                params.field = Field.SINE;
                // modify parameters here for interactive testing
            }

            MappedSingletDensity sim = new MappedSingletDensity(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            int blockSize = 100;
            MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
            densityMeter.setBox(sim.box());
            densityMeter.setProfileDim(2);
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            densityMeter.setDataSource(meterNMolecules);
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(5 * blockSize);
            DataPumpListener pump = new DataPumpListener(densityMeter, acc, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pump);

            MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
            densityMeterForce.setProfileDim(2);
            AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpForce);

            FunctionDifferentiable f;
            double L = sim.box().getBoundary().getBoxSize().getX(2);
            switch (params.field) {
                case SINE:
                    f = new FunctionSine(5, L);
                    break;
                case UNIFORM:
                    f = new FunctionUniform(L);
                    break;
                default:
                    throw new RuntimeException("not yet");
            }
            MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterMappedAvg.setProfileDim(2);
            AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);

            MeterProfileMappedAvg densityMeterP = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterP.getXDataSource().setNValues(104);
            densityMeterP.setProfileDim(2);
            densityMeterP.setBehavior(MeterProfileMappedAvg.Behavior.P);
            MeterProfileMappedAvg densityMeterZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot0.getXDataSource().setNValues(104);
            densityMeterZidot0.setProfileDim(2);
            densityMeterZidot0.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot1.getXDataSource().setNValues(104);
            densityMeterZidot1.setProfileDim(2);
            densityMeterZidot1.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot2.getXDataSource().setNValues(104);
            densityMeterZidot2.setProfileDim(2);
            densityMeterZidot2.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot2.setZidotZ(L / 4);

            MeterProfileMappedAvg densityMeterDZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot0.getXDataSource().setNValues(104);
            densityMeterDZidot0.setProfileDim(2);
            densityMeterDZidot0.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterDZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot1.getXDataSource().setNValues(104);
            densityMeterDZidot1.setProfileDim(2);
            densityMeterDZidot1.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterDZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot2.getXDataSource().setNValues(104);
            densityMeterDZidot2.setProfileDim(2);
            densityMeterDZidot2.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot2.setZidotZ(L / 4);

            DisplayPlot densityPlot = new DisplayPlot();
            acc.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accForce.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accMappedAvg.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});

            new DataPump(densityMeterP, densityPlot.getDataSet().makeDataSink()).actionPerformed();

            densityPlot.setLabel("density");
            graphic.add(densityPlot);

            DisplayPlot errorPlot = new DisplayPlot();
            acc.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accForce.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accMappedAvg.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            errorPlot.setLabel("error");
            graphic.add(errorPlot);

            DisplayPlot corPlot = new DisplayPlot();
            acc.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accForce.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accMappedAvg.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            corPlot.setLabel("correlation");
            graphic.add(corPlot);

            DisplayPlot zidotPlot = new DisplayPlot();
            new DataPump(densityMeterZidot0, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot1, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot2, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            zidotPlot.setLabel("zidot");
            graphic.add(zidotPlot);

            DisplayPlot dzidotPlot = new DisplayPlot();
            new DataPump(densityMeterDZidot0, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot1, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot2, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            dzidotPlot.setLabel("-dzidot");
            graphic.add(dzidotPlot);

            List<DataPump> pumps = graphic.getController().getDataStreamPumps();
            pumps.add(pump);
            pumps.add(pumpForce);
            pumps.add(pumpMappedAvg);

            graphic.makeAndDisplayFrame();

        }
    }

    enum Field {
        SINE, UNIFORM, PARABOLIC
    }

    public static class SimParams extends ParameterBase {
        public long steps = 1000000;
        public double density = 0.35;
        public double temperature = 3;
        public int numAtoms = 500;
        public Field field = Field.UNIFORM;
    }
}
