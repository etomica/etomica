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
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestLJMD3D extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;

    public TestLJMD3D(int numAtoms, int numSteps, Configuration config) {
        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, 4, space);
        double sigma = 1.0;
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, potential, 3);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(p2, new AtomType[]{leafType, leafType});

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3D.class);

        TestLJMD3D sim = new TestLJMD3D(numAtoms, params.numSteps / params.numAtoms, config);

        MeterPressure pMeter = new MeterPressure(sim.getSpace());
        pMeter.setIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(10);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 4);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pPump);
        pumpListener.setInterval(4);
        sim.integrator.getEventManager().addListener(pumpListener);
        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
        energyMeter.setBox(sim.box);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPumpListener energyManager = new DataPumpListener(energyMeter, energyAccumulator, 4);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyManager));

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.currentTimeMillis();

        double Z = ((DataDouble) ((DataGroup) pAccumulator.getData()).getData(pAccumulator.AVERAGE.index)).x * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x;
        avgPE /= numAtoms;
        System.out.println("Z=" + Z);
        System.out.println("PE/epsilon=" + avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k=" + Cv);

        if (Double.isNaN(Z) || Math.abs(Z + 0.25) > 0.2) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE + 4.56) > 0.06) {
            System.exit(1);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv - 0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(1);
        }
        System.out.println("runtime: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 1000000;
    }

}
