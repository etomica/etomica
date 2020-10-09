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
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureFromIntegratorFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestLJMD3D extends Simulation {

    public IntegratorVelocityVerletFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;

    public TestLJMD3D(int numAtoms, Configuration config) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox();
        PotentialMasterListFasterer potentialMaster = new PotentialMasterListFasterer(this, box, 2, 4);
        double sigma = 1.0;
        integrator = new IntegratorVelocityVerletFasterer(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, potential, 3);
        AtomType leafType = species.getLeafType();

        potentialMaster.setPairPotential(leafType, leafType, p2);

        integrator.getEventManager().addListener(potentialMaster);

        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", numAtoms), TestLJMC3DSlowerer.class);

        TestLJMD3D sim = new TestLJMD3D(numAtoms, config);

        MeterPressureFromIntegratorFasterer pMeter = new MeterPressureFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(50);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 4);
        sim.integrator.getEventManager().addListener(pPump);
        MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(50);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 4);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / params.numAtoms));
        long t2 = System.currentTimeMillis();

        double Z = ((DataDouble) ((DataGroup) pAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        avgPE /= numAtoms;
        System.out.println("Z=" + Z);
        System.out.println("PE/epsilon=" + avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k=" + Cv);

        if (Double.isNaN(Z) || Math.abs(Z + 0.25 - 0.4) > 0.2) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE + 4.56 - 0.2) > 0.06) {
            System.exit(2);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv - 0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(3);
        }
        System.out.println("runtime: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 1000000;
    }

}
