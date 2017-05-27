/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.graphics.SimulationGraphic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;

/**
 * Simulation for rheology module.
 * 
 * @author Andrew Schultz
 */
public class SimulationRheology extends Simulation {

    public final Box box;
    public final SpeciesSpheres species;
    public final IntegratorPolymer integrator;
    public final ActivityIntegrate activityIntegrate;
    public final ConformationPolymer conformation;
    
    public SimulationRheology(Space space) {
        super(space);
        box = new Box(new BoundaryRectangularNonperiodic(space), space);
        Vector d = space.makeVector();
        d.E(20);
        box.getBoundary().setBoxSize(d);
        addBox(box);
        species = new SpeciesSpheres(this, space, 2);
        species.setIsDynamic(true);
        addSpecies(species);
        box.setNMolecules(species, 1);
        conformation = new ConformationPolymer(space, random);
        conformation.initializePositions(box.getMoleculeList().getMolecule(0).getChildList());
        integrator = new IntegratorPolymer(null, getRandom(), 0.01, 1.0, space);
        integrator.setB(1);
        integrator.setBox(box);
        activityIntegrate = new ActivityIntegrate(integrator, 0, false);
        getController().addAction(activityIntegrate);
    }
    
    public void setChainLength(int newChainLength) {
        if (newChainLength < 2) {
            throw new IllegalArgumentException("too short");
        }
        box.setNMolecules(species, 0);
        species.setNumLeafAtoms(newChainLength);
        box.setNMolecules(species, 1);
        conformation.initializePositions(box.getMoleculeList().getMolecule(0).getChildList());
    }

    public int getChainLength() {
        return species.getNumLeafAtoms();
    }
    
    private static final long serialVersionUID = 1L;

    public static void main(String[] args) {
        SimulationRheology sim = new SimulationRheology(Space3D.getInstance());
        sim.setChainLength(10);
        SimulationGraphic graphic = new SimulationGraphic(sim, sim.getSpace(), sim.getController());
        graphic.setPaintInterval(sim.box, 1);
        graphic.makeAndDisplayFrame();
    }
}
