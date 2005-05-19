/*
 * Created on Dec 3, 2004
 */
package etomica.models.hexane;

import etomica.ConfigurationLattice;
import etomica.Space;


/**
 * @author nancycribbin
 *
 * Creates the lattice that the hexane molecules of Dr. Monson's data are placed
 * into.
 */
public class ConfigurationHexane extends ConfigurationLattice {

	ConfigurationHexane(Space space) {
	    super(new PrimitiveHexane(space));
	}
	
//	//DUDE,  we need to figure out how to write this... 
//	public void initializePositions(AtomList[] lists){
//	    //Do nothing if there are no atoms on the list.
//	    if(lists.length == 0) {return;}
//	    
//	    atomIterator.setLists(lists);
//	    
//	    //Sequential has dimension stuff here.
//	    
//	    int sumOfMolecules = atomIterator.size();
//	    if(sumOfMolecules == 0) {return;}
//	    
//	    //We know we are in 3 dimensions, so we don't have to worry about putting in a lineLattice the way sequential does
//	    
//	    
//	    
//	    //Place molecules
//	    int i = 0;
//	    atomIterator.reset();
//	    while(atomIterator.hasNext()){
//	        Atom a = atomIterator.nextAtom();
//	        //initialize coordingats of child atoms
//	        if(!a.node.isLeaf()){
//	            Conformation conform = a.type.creator().getConformation();
//	            conform.initializePositions(((AtomTreeNodeGroup)a.node).childlist);
//	            
//	        } else {
//	            a.coord.position().E(rLat[i]);
//	        }
//	        i++
//	    }
//	    
//	}
	
	
}