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
import etomica.potential.IPotentialAtomic;
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
    public final IPotentialAtomic p1;

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
            p1 = new P1Sine(space, 5, params.temperature);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        }
        else {
            p1 = null;
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
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pump = new DataPumpListener(densityMeter, acc, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pump);

            MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
            densityMeterForce.setProfileDim(2);
            AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpForce);

            FunctionDifferentiable f;
            switch (params.field) {
                case SINE:
                    f = new FunctionSine(5, sim.box().getBoundary().getBoxSize().getX(2));
                    break;
                default:
                    throw new RuntimeException("not yet");
            }
            MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterMappedAvg.setProfileDim(2);
            AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);


            DisplayPlot densityPlot = new DisplayPlot();
            acc.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accForce.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accMappedAvg.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            densityPlot.setLabel("density");
            graphic.add(densityPlot);

            DisplayPlot errorPlot = new DisplayPlot();
            acc.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accForce.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accMappedAvg.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            errorPlot.setLabel("error");
            graphic.add(errorPlot);

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
