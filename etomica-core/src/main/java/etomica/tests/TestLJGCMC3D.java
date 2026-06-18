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
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
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
 * Grand Canonical Lennard-Jones Monte Carlo simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 * This class uses the same LJMC3D that TestLJMC3DSlowerer uses
 */
public class TestLJGCMC3D extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public MCMoveInsertDelete mcMoveID;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;

    public TestLJGCMC3D(int numAtoms, Configuration config) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        box = this.makeBox();

        PotentialMasterCell potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 3, BondingInfo.noBonding());
        potentialMaster.doOneTruncationCorrection = true;
        double sigma = 1.0;
        integrator = new IntegratorMC(potentialMaster, this.getRandom(), 1.0, box);
        integrator.setTemperature(1.1);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        mcMoveAtom.setStepSize(0.2 * sigma);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setTunable(false);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveID = new MCMoveInsertDelete(potentialMaster, random, space);
        mcMoveID.setBox(box);
        mcMoveID.setMu(-3.67);
        integrator.getMoveManager().addMCMove(mcMoveID);
        integrator.getMoveManager().setEquilibrating(false);
        mcMoveID.setSpecies(species);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(sigma, 1.0);
        double truncationRadius = 3.0 * sigma;
        if (truncationRadius > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
        AtomType leafType = species.getLeafType();
        potentialMaster.setPairPotential(leafType, leafType, potentialTruncated);

        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJGCMC3D.class);

        TestLJGCMC3D sim = new TestLJGCMC3D(numAtoms, config);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));

        int pInterval = 2 * numAtoms;
        int bs = params.numSteps / (pInterval * 50);
        if (bs == 0) bs = 1;
        MeterPressure pMeter = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
        pMeter.setTemperature(sim.integrator.getTemperature());
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, pInterval);
        sim.integrator.getEventManager().addListener(pPump);

        bs = params.numSteps / 50;
        if (bs == 0) bs = 1;
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener uPump = new DataPumpListener(energyMeter, energyAccumulator);
        sim.integrator.getEventManager().addListener(uPump);

        MeterDensity densityMeter = new MeterDensity(sim.box);
        AccumulatorAverage densityAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pumpDensity = new DataPumpListener(densityMeter, densityAccumulator);
        sim.integrator.getEventManager().addListener(pumpDensity);

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.currentTimeMillis();

        System.out.println("runtime: " + (t2 - t1) * 0.001);

        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + avgP + " " + errP + " " + corP);

        double avgRho = densityAccumulator.getData(densityAccumulator.AVERAGE).getValue(0);
        double errRho = densityAccumulator.getData(densityAccumulator.ERROR).getValue(0);
        double corRho = densityAccumulator.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
        double numAtomsAvg = avgRho * sim.box.getBoundary().volume();

        double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtomsAvg;
        double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtomsAvg;
        double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("PE " + avgPE + " " + errPE + " " + corPE);

        System.out.println("rho " + avgRho + " " + errRho + " " + corRho);

        // expected values based on 10^8 steps
        // stdev based on 10^8 steps for N=500 uncertainty, scaled up to 10^6 steps
        //   other sims have significant correlation
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.0815; // finite size effect smaller than uncertainty
        double stdevP = 0.03;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedPE = -4.492; // finite size effect comparable to uncertainty
        double stdevPE = 0.09;
        if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }

        double expectedRho = 0.65; // finite size effect smaller than uncertainty
        double stdevRho = 0.006;
        if (Double.isNaN(avgRho) || Math.abs(avgRho - expectedRho) / stdevRho > 4) {
            System.exit(3);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 1000000;
    }
}
