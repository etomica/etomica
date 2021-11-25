/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Iron;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureFromIntegrator;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.BondingInfo;
import etomica.potential.EmbeddingSqrt;
import etomica.potential.P2SoftSphereFloatTab;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.compute.PotentialComputeEAM;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * EAM Iron simulation with an FCC crystal.
 */
public class TestEAM extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesGeneral species;
    public Box box;

    public TestEAM(int numAtoms, double temperature) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.element(Iron.INSTANCE), true);
        addSpecies(species);

        box = this.makeBox();
        NeighborListManager nbrs = new NeighborListManager(this.getSpeciesManager(), box, 2, 7.2, BondingInfo.noBonding());
        nbrs.setDoDownNeighbors(true);
        PotentialComputeEAM potentialMaster = new PotentialComputeEAM(getSpeciesManager(), box, nbrs);
        potentialMaster.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.001, temperature, box);
        integrator.setIsothermal(true);

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.15);
        inflater.actionPerformed();

        AtomType leafType = species.getLeafType();
        P2SoftSphereFloatTab p2 = new P2SoftSphereFloatTab(1, 1.2446188036318708E7, 8.7932, 6, 1000);
        P2SoftSphericalTruncated p2t = new P2SoftSphericalTruncated(p2, 6);
        potentialMaster.setPairPotential(leafType, leafType, p2t);
        P2SoftSphereFloatTab pRho = new P2SoftSphereFloatTab(1, 26068.513192447575, 8.14475, 6, 1000);
        P2SoftSphericalTruncated pRhot = new P2SoftSphericalTruncated(pRho, 6);
        potentialMaster.setRhoPotential(leafType, pRhot);
        EmbeddingSqrt f = new EmbeddingSqrt(6129.634374295454);
        potentialMaster.setEmbeddingPotential(leafType, f);

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;

        TestEAM sim = new TestEAM(numAtoms, Kelvin.UNIT.toSim(params.temperatureK));
        sim.integrator.reset();
        System.out.println("u0: " + sim.integrator.getPotentialEnergy() / numAtoms);

        int steps = params.numSteps / numAtoms;
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));

        int bs = steps / (4 * 100);
        MeterPressureFromIntegrator pMeter = new MeterPressureFromIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 4);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 4);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.currentTimeMillis();

        System.out.println("runtime: " + (t2 - t1) * 0.001);

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
        // for MD, avg values are very close for short and longer runs
        // stdev based on 50 x 10^6 steps with 4000 atoms (a bit larger than for 500)
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.466751 + 9.74 / numAtoms;
        double stdevP = 0.005;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedPE = -3.81607 + 4.40 / numAtoms;
        double stdevPE = 0.0012;
        if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }

        double expectedCv = 0.3946 + 7.724 / numAtoms;
        double stdevCv = 0.038; // stdev 500 atoms is ~2x smaller
        // at 4sigma, this isn't too useful expect that it's not super-big
        if (Double.isNaN(Cv) || Math.abs(Cv - expectedCv) / stdevCv > 4) {
            System.exit(3);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 864;
        public int numSteps = 10_000_000;
        public double temperatureK = 6000;
    }

}
