/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
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

public class TestHSMD3DSlow extends Simulation {

    public IntegratorHard integrator;
    public SpeciesGeneral species, species2;
    public Box box;

    public TestHSMD3DSlow(Space _space, int numAtoms, Configuration config) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species2);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);
        box = makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type1, type1});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type1, type2});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type2, type2});
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
        Configuration config = Configurations.fromResourceFile(String.format("HSMD3D%d.pos", numAtoms), TestHSMD3DSlow.class);

        TestHSMD3DSlow sim = new TestHSMD3DSlow(Space3D.getInstance(), numAtoms, config);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / numAtoms / 10));
        sim.integrator.resetStepCount();
        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);

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
        public int numSteps = 20000000;
    }
}
