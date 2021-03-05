/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 *
 * @author David Kofke
 */

public class TestSWMD3DSlow extends Simulation {

    public IntegratorHard integrator;
    public SpeciesGeneral species, species2;
    public Box box;

    public TestSWMD3DSlow(Space _space, int numAtoms, Configuration config) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species2);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        double neighborRangeFac = 1.9;
        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);
        box = makeBox();
        integrator = new IntegratorHard(random, potentialMaster, box);
        integrator.setTimeStep(0.005);
        integrator.setIsothermal(true);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();

        potentialMaster.addPotential(new P2SquareWell(space, sigma, 1.5, 0.5, false), new AtomType[]{type1, type1});

        potentialMaster.addPotential(new P2SquareWell(space, sigma, 1.5, 0.5, false), new AtomType[]{type1, type2});

        potentialMaster.addPotential(new P2SquareWell(space, sigma, 1.5, 0.5, false), new AtomType[]{type2, type2});
        box.setNMolecules(species, numAtoms);
        box.setNMolecules(species2, numAtoms / 100);
        box.getBoundary().setBoxSize(Vector.of(l, l, l));
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        config.initializeCoordinates(box);

//        WriteConfiguration writeConfig = new WriteConfiguration("foo",box,1);
//        integrator.addIntervalListener(writeConfig);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("HSMD3D%d.pos", numAtoms), TestSWMD3DSlow.class);

        TestSWMD3DSlow sim = new TestSWMD3DSlow(Space3D.getInstance(), numAtoms, config);
        long steps = params.numSteps / numAtoms;
//        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
//        sim.integrator.resetStepCount();

        int interval = 4;
        long bs = steps / (interval * 100);
        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, interval);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegrator uMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage uAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener uPump = new DataPumpListener(uMeter, uAccumulator, interval);
        sim.integrator.getEventManager().addListener(uPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);

        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + avgP + " " + errP + " " + corP);

        double avgU = uAccumulator.getData(uAccumulator.AVERAGE).getValue(0) / numAtoms;
        double errU = uAccumulator.getData(uAccumulator.ERROR).getValue(0) / numAtoms;
        double corU = uAccumulator.getData(uAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("U " + avgU + " " + errU + " " + corU);

        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) uAccumulator.getData()).getData(uAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k " + Cv);

        // expected values based on 2x10^9 / numAtoms steps for N=500 and 4000 (32000 not tested)
        // equilibration (initial config is HS) and block size are minor issues for N=4000
        // for MD, avg values are very close for short and longer runs
        // stdev based on 100 x (2x10^7/numAtoms) steps with 500 atoms (a bit smaller for P with 4000 atoms)
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 1.807199 - 0.807742 / numAtoms;
        double stdevP = 0.008;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedU = -2.497658 + 0.128618 / numAtoms;
        double stdevU = 0.0009;
        if (Double.isNaN(avgU) || Math.abs(avgU - expectedU) / stdevU > 4) {
            System.exit(2);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 50000000;
    }
}
