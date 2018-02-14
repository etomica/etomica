/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Box box;

    public TestHSMD3D(Space _space, int numAtoms, int numSteps, Configuration config) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573 * Math.pow((numAtoms / 2000.0), 1.0 / 3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);
        box = new Box(space);
        integrator = new IntegratorHard(this, potentialMaster, space, box);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(numSteps);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        removeSpecies(species);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        species2 = new SpeciesSpheresMono(this, space);
        species2.setIsDynamic(true);
        addSpecies(species2);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type1, type1});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type1, type2});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false), new AtomType[]{type2, type2});
        addBox(box);
        box.setNMolecules(species, numAtoms);
        box.setNMolecules(species2, numAtoms / 100);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{l, l, l}));
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
        ConfigurationFile config = new ConfigurationFile("HSMD3D"+Integer.toString(numAtoms));

        TestHSMD3D sim = new TestHSMD3D(Space3D.getInstance(), numAtoms, params.numSteps / numAtoms, config);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);

        sim.getController().actionPerformed();

        double Z = pMeter.getDataAsScalar()*sim.box.getBoundary().volume()/(sim.box.getMoleculeList().getMoleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);

        // compressibility factor for this system should be 5.22
        if (Double.isNaN(Z) || Math.abs(Z-5.22) > 0.04) {
            System.exit(1);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 20000000;
    }
}
