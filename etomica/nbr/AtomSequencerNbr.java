/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomArrayList;
import etomica.AtomLinker;
import etomica.AtomSequencerFactory;
import etomica.Potential;
import etomica.nbr.cell.AtomSequencerCell;
import etomica.utility.Arrays;

/**
 * Sequencer used to maintain neighbor lists.  Holds lists of atoms
 * that were elsewhere deemed to be neighbors of the sequencer's atom.
 * Atoms are stored with reference to the potential that governs their
 * interactions.
 */
public class AtomSequencerNbr extends AtomSequencerCell {

	protected AtomArrayList[] upList, downList;
	
	/**
	 * Constructs sequencer for the given atom.
	 */
	public AtomSequencerNbr(Atom a) {
		super(a);
		upList = new AtomArrayList[0];
		downList = new AtomArrayList[0];
	}
	
    /**
     * Adds the given atom to the neighbor set that are uplist of
     * this sequencer's atom and which interact via the given potential.  
     * @param a the new uplist neighbor atom
     * @param potential the potential between the atoms
     */
	public void addUpNbr(Atom a, Potential potential) {
		int index = 0;
		try {
			index = atom.type.getNbrManagerAgent().getPotentialIndex(potential);
			upList[index].add(a);
		} catch(ArrayIndexOutOfBoundsException e) {
			index = addPotential(potential);
			upList[index].add(a);
		}
	}

    /**
     * Adds the given atom to the neighbor set that are downlist of
     * this sequencer's atom and which interact via the given potential.  
     * @param a the new downlist neighbor atom
     * @param potential the potential between the atoms
     */
    public void addDownNbr(Atom a, Potential potential) {
		int index = 0;
		try {
			index = atom.type.getNbrManagerAgent().getPotentialIndex(potential);
			downList[index].add(a);
		} catch(ArrayIndexOutOfBoundsException e) {
			index = addPotential(potential);
			downList[index].add(a);
		}
	}
	
    /**
     * Returns an array of uplist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential.
     */
	public AtomArrayList[] getUpList() {
		return upList;
	}
	
    /**
     * Returns an array of downlist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential.
     */
	public AtomArrayList[] getDownList() {
		return downList;
	}
	
    /**
     * Indicates that the given potential is involved in this sequencer's
     * atom's interactions.  Enables lists to be stored of atoms interacting 
     * with this one via the given potential.
     * @param p
     * @return index of the neighbor-list arrays giving the list of atoms
     * associated with the potential
     */
	public int addPotential(Potential p) {
		int index = atom.type.getNbrManagerAgent().addPotential(p);
		if (index > upList.length-1) {
            upList = (AtomArrayList[])Arrays.addObject(upList, new AtomArrayList());
            downList = (AtomArrayList[])Arrays.addObject(downList, new AtomArrayList());
		}
		return index;
	}
	
    /**
     * Indicates that neighbor lists will no longer be kept for the given potential.
     * @param p
     */
    //TODO consider whether this method might foul up the index assignments
	public void removePotential(Potential p) {
		atom.type.getNbrManagerAgent().removePotential(p);
		if (upList.length == 0) throw new RuntimeException("potential list empty in removePotential");
		upList = new AtomArrayList[upList.length-1];
		downList = new AtomArrayList[downList.length-1];
		for (int i=0; i<upList.length; i++) {
			upList[i] = new AtomArrayList();
			downList[i] = new AtomArrayList();
		}
	}
	
    /**
     * Clears neighbor lists, removing all listed neighbor atoms.
     */
	public void clearNbrs() {
		int length = upList.length;
		for (int i=0; i<length; i++) {
			upList[i].clear();
			downList[i].clear();
		}
	}

    /**
     * A factory class that will make this atom sequencer.  Typically this
     * is passed to the constructor of the atom factory.
     */
    public static final AtomSequencerFactory FACTORY = new AtomSequencerFactory() {
        public AtomLinker makeSequencer(Atom atom) {
            return new AtomSequencerNbr(atom);
        }
    };
    
}
