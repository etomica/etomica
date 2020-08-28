/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;


import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation for rheology module.
 * 
 * @author Andrew Schultz
 */
public class SimulationRheology extends Simulation {

    public final Box box;
    public final SpeciesSpheresMono species;
    public final IntegratorPolymer integrator;

    public final ConfigurationPolymer config;
    
    public SimulationRheology(Space space) {
        super(space);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        Vector d = space.makeVector();
        d.E(20);
        box.getBoundary().setBoxSize(d);
        box.setNMolecules(species, 2);
        config = new ConfigurationPolymer(space, random);
        config.initializeCoordinates(box);
        integrator = new IntegratorPolymer(null, getRandom(), 0.01, 1.0, box);
        integrator.setB(1);
        getController().addActivity(new ActivityIntegrate(integrator, true));
    }
    
    public void setChainLength(int newChainLength) {
        if (newChainLength < 2) {
            throw new IllegalArgumentException("too short");
        }
        box.setNMolecules(species, newChainLength);
        config.initializeCoordinates(box);
    }

    public int getChainLength() {
        return box.getNMolecules(species);
    }

    public static void main(String[] args) {
        SimulationRheology sim = new SimulationRheology(Space3D.getInstance());
        sim.setChainLength(10);
        SimulationGraphic graphic = new SimulationGraphic(sim);
        graphic.setPaintInterval(sim.box, 1);
        graphic.makeAndDisplayFrame();
    }
}
