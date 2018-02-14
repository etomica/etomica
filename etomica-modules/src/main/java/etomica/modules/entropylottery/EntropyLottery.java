/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import etomica.action.activity.ActivityIntegrate;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.species.SpeciesSpheresMono;

public class EntropyLottery extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    
    public EntropyLottery(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster();
        
        final int N = 30;  //number of atoms
        
        //controller and integrator
        integrator = new IntegratorMC(this, potentialMaster, box);
        MCMoveAtomAdjacent move = new MCMoveAtomAdjacent(getRandom(), space);
        integrator.getMoveManager().addMCMove(move);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
        //construct box
	    box = new Box(new BoundaryRectangularNonperiodic(space), space);
        addBox(box);
        box.setNMolecules(species, N);
        Vector dimensions = space.makeVector();
        dimensions.E(10);
        box.getBoundary().setBoxSize(dimensions);
        new ConfigurationZero(space).initializeCoordinates(box);
        integrator.setBox(box);
		
//        BoxImposePbc imposePBC = new BoxImposePbc(box);
//        integrator.addListener(new IntervalActionAdapter(imposePBC));
        
    }  
    
    public static void main(String[] args) {
        EntropyLottery sim = new EntropyLottery(Space1D.getInstance());
        sim.getController().actionPerformed();
    }
    
}


