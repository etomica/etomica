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
import etomica.data.meter.MeterPressureHardFasterer;
import etomica.integrator.IntegratorHardFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.PotentialComputePair;
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
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {

    public IntegratorHardFasterer integrator;
    public SpeciesGeneral species, species2;
    public Box box;

    public TestHSMD3D(Space _space, int numAtoms, Configuration config) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species2);

        box = makeBox();

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 1, 1.5, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();

        potentialMaster.setPairPotential(type1, type1, P2HardSphere.makePotential(sigma));
        potentialMaster.setPairPotential(type1, type2, P2HardSphere.makePotential(sigma));
        potentialMaster.setPairPotential(type2, type2, P2HardSphere.makePotential(sigma));

        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager, random, 0.01, 1.0, box);

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
        Configuration config = Configurations.fromResourceFile(String.format("HSMD3D%d.pos", numAtoms), TestHSMD3D.class);

        TestHSMD3D sim = new TestHSMD3D(Space3D.getInstance(), numAtoms, config);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / numAtoms / 10));
        sim.integrator.resetStepCount();

        MeterPressureHardFasterer pMeter = new MeterPressureHardFasterer(sim.integrator);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / numAtoms));
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);

        double Z = pMeter.getDataAsScalar() * sim.box.getBoundary().volume() / (sim.box.getMoleculeList().size() * sim.integrator.getTemperature());
        System.out.println("Z=" + Z);

        // compressibility factor for this system should be 5.22
        if (Double.isNaN(Z) || Math.abs(Z - 5.22) > 0.04) {
            System.exit(1);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 50000000;
    }
}
