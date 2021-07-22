/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureHardFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHardFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2SquareWell;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple square-well molecular dynamics simulation in 3D.
 * This uses initial configs from TestHSMD3D, but runs SQW at high enough T so
 * that lack of equilibration is not too bad.
 */

public class TestSWMD3D extends Simulation {

    public IntegratorHardFasterer integrator;
    public SpeciesGeneral species, species2;
    public Box box;

    public TestSWMD3D(Space _space, int numAtoms, Configuration config) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species2);

        box = makeBox();

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 1, 1.8, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();

        potentialMaster.setPairPotential(type1, type1, P2SquareWell.makePotential(sigma, 1.5, 0.5));
        potentialMaster.setPairPotential(type1, type2, P2SquareWell.makePotential(sigma, 1.5, 0.5));
        potentialMaster.setPairPotential(type2, type2, P2SquareWell.makePotential(sigma, 1.5, 0.5));
        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager, random, 0.005, 1.0, box, getSpeciesManager());
        integrator.setIsothermal(true);
        box.setNMolecules(species, numAtoms);
        box.setNMolecules(species2, numAtoms / 100);
        box.getBoundary().setBoxSize(Vector.of(l, l, l));
        config.initializeCoordinates(box);
        new BoxImposePbc(box, space).actionPerformed();
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        long steps = params.numSteps / numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("HSMD3D%d.pos", numAtoms), TestSWMD3D.class);

        TestSWMD3D sim = new TestSWMD3D(Space3D.getInstance(), numAtoms, config);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        sim.integrator.resetStepCount();

        int interval = 4;
        long bs = steps / (interval * 100);
        MeterPressureHardFasterer pMeter = new MeterPressureHardFasterer(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, interval);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegratorFasterer uMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
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

        double avgU = uAccumulator.getData(pAccumulator.AVERAGE).getValue(0) / numAtoms;
        double errU = uAccumulator.getData(pAccumulator.ERROR).getValue(0) / numAtoms;
        double corU = uAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("U " + avgU + " " + errU + " " + corU);

        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) uAccumulator.getData()).getData(uAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k " + Cv);

        // expected values based on 2x10^9 / numAtoms steps for N=500 and 4000 (32000 not tested)
        // equilibration (initial config is HS) and block size are minor issues for N=4000
        // for MD, avg values are very close for short and longer runs
        // stdev based on 100 x (2x10^7/numAtoms) steps (max of N=500 and 4000 data)
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 1.807199 - 0.807742 / numAtoms;
        double stdevP = 0.004;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedU = -2.497658 + 0.128618 / numAtoms;
        double stdevU = 0.0006;
        if (Double.isNaN(avgU) || Math.abs(avgU - expectedU) / stdevU > 4) {
            System.exit(2);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 50000000;
    }
}
