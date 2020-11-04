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
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureFromIntegratorFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.ewald.P2Ewald6Real;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestLJMD3DEwald extends Simulation {

    public IntegratorVelocityVerletFasterer integrator;
    public SpeciesGeneral species;
    public Box box;

    public TestLJMD3DEwald(int numAtoms, Configuration config) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox();
        PotentialMasterListFasterer pair = new PotentialMasterListFasterer(this, box, 2, 4, BondingInfo.noBonding());
        PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(this, box, BondingInfo.noBonding());
        PotentialComputeAggregate aggregate = new PotentialComputeAggregate(pair, ewaldFourier);
//        PotentialComputeAggregate aggregate = new PotentialComputeAggregate(pair);
        integrator = new IntegratorVelocityVerletFasterer(aggregate, random, 0.01, 1.1, box);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        double alpha6 = 1;
        P2Ewald6Real ewaldReal = new P2Ewald6Real(1, 1, 1, 1, 1);
        P2SoftSphere pCore12 = new P2SoftSphere(space, 1, 4, 12);
        P2SoftSphericalTruncatedSum2 trunc = new P2SoftSphericalTruncatedSum2(space, ewaldReal, pCore12, 5);
        AtomType leafType = species.getLeafType();

        pair.setPairPotential(leafType, leafType, trunc);
        ewaldFourier.setAlpha6(alpha6);
        ewaldFourier.setkCut(1.5);
        ewaldFourier.setAlpha(0);
        ewaldFourier.setR6Coefficient(leafType, 1, 1);

        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3DSlowerer.class);

        TestLJMD3DEwald sim = new TestLJMD3DEwald(numAtoms, config);

        sim.integrator.reset();
        double u = sim.integrator.getPotentialCompute().computeAll(false);
        System.out.println(u / numAtoms);
        System.out.println(sim.integrator.getPotentialCompute().getLastVirial());

        int steps = params.numSteps / params.numAtoms;
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));

        int bs = steps / (4 * 100);
        MeterPressureFromIntegratorFasterer pMeter = new MeterPressureFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 4);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 4);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / params.numAtoms));
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
        public int numAtoms = 500;
        public int numSteps = 5000000;
    }

}
