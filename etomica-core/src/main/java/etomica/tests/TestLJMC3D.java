/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class TestLJMC3D extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public PotentialMasterCell potentialMaster;
    public P2LennardJones potential;

    public TestLJMC3D(int numAtoms, Configuration config) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        double sigma = 1.0;
        box = this.makeBox();
        potentialMaster = new PotentialMasterCell(this.getSpeciesManager(), box, 2, BondingInfo.noBonding());
        integrator = new IntegratorMC(potentialMaster, random, 1.1, box);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        mcMoveAtom.setStepSize(0.275 * sigma);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setTunable(false);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().setEquilibrating(false);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);

        double truncationRadius = 3.0 * sigma;

        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        AtomType leafType = species.getLeafType();
        potentialMaster.setPairPotential(leafType, leafType, potentialTruncated);

        config.initializeCoordinates(box);

        potentialMaster.init();
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3D.class);

        TestLJMC3D sim = new TestLJMC3D(numAtoms, config);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));

        int bs = params.numSteps / (100 * 2 * numAtoms);
        MeterPressure pMeter = new MeterPressure(sim.box, sim.potentialMaster);
        pMeter.setTemperature(sim.integrator.getTemperature());
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 2 * numAtoms);
        sim.integrator.getEventManager().addListener(pPump);

        bs = params.numSteps / (100 * 10);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 10);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);

        System.out.println("Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());

        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + avgP + " " + errP + " " + corP);

        double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtoms;
        double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtoms;
        double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("PE " + avgPE + " " + errPE + " " + corPE);

        double temp = sim.integrator.getTemperature();
        double Cv = energyAccumulator.getData(energyAccumulator.STANDARD_DEVIATION).getValue(0);
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv " + Cv);

        // expected values based on 10^8 steps
        // stdev based on 50 x 10^6 steps with 4000 atoms (a bit larger than for 500)
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.0826 - 3.71 / numAtoms;
        double stdevP = 0.015;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedPE = -4.48971 - 0.331 / numAtoms;
        double stdevPE = 0.004;
        if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }

        double expectedCv = 0.597;
        double stdevCv = 0.1; // stdev 500 atoms is ~2x smaller
        // at 4sigma, this isn't too useful expect that it's not super-big
        if (Double.isNaN(Cv) || Math.abs(Cv - expectedCv) / stdevCv > 4) {
            System.exit(3);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 1000000;
    }
}
