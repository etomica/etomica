/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Mesoscale simulation for Droplet module.
 *
 * @author Andrew Schultz
 */
public class Droplet extends Simulation {

    public final SpeciesSpheresMono species;
    public final Box box;
    public final IntegratorDroplet integrator;
    public final ActivityIntegrate activityIntegrate;
    public final P2Cohesion p2;
    public final P1Smash p1Smash;
    public final ConfigurationDroplet config;
    public final AtomFilterLiquid liquidFilter;
    public final MeterDeformation meterDeformation;

    public Droplet(Space _space) {
        super(_space);

        //species
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        int numAtoms = 2000;
        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this);

        //controller and integrator
        integrator = new IntegratorDroplet(this, potentialMaster, box);
        activityIntegrate = new ActivityIntegrate(integrator);
//        activityIntegrate.setMaxSteps(10);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.2);
        integrator.setTemperature(0);

        //potentials
        AtomType leafType = species.getLeafType();

        p2 = new P2Cohesion(space);
        p2.setEpsilon(1.0);
        double vol = 4.0 / 3.0 * Math.PI;
        p2.setDv(vol / numAtoms);
        potentialMaster.addPotential(p2, new AtomType[]{leafType, leafType});

        p1Smash = new P1Smash(space);
        p1Smash.setG(1.0);
        potentialMaster.addPotential(p1Smash, new AtomType[]{leafType});

        //construct box
        Vector dim = space.makeVector();
        dim.E(new double[]{4, 4, 4});
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, numAtoms);

        config = new ConfigurationDroplet(random, space);
        config.initializeCoordinates(box);

        meterDeformation = new MeterDeformation(space);
        meterDeformation.setBox(box);

        liquidFilter = new AtomFilterLiquid(space, meterDeformation);
        liquidFilter.setCutoff(0.9);

        p2.setLiquidFilter(liquidFilter);
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        Droplet sim = new Droplet(space);
        sim.getController().actionPerformed();
    }//end of main
}
