/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbratom;

import etomica.Atom;
import etomica.Potential;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;
import etomica.nbratom.cell.AtomSequencerCell;
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
        int index = atom.type.getNbrManagerAgent().getPotentialIndex(potential);
        ensureCapacity(index);
        upList[index].add(a);
//        if (atom.node.getOrdinal() == 1) {
//            System.out.println("in aUN, "+index+" "+upList[index].size());
//        }
    }

    /**
     * Adds the given atom to the neighbor set that are downlist of
     * this sequencer's atom and which interact via the given potential.  
     * @param a the new downlist neighbor atom
     * @param potential the potential between the atoms
     */
    public void addDownNbr(Atom a, Potential potential) {
        int index = atom.type.getNbrManagerAgent().getPotentialIndex(potential);
        ensureCapacity(index);
        downList[index].add(a);
//        if (atom.node.getOrdinal() == 1) {
//            System.out.println("in aDN, "+index+" "+downList[index].size());
//        }
    }

    /**
     * Returns an array of uplist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
     */
    public AtomArrayList[] getUpList() {
        return upList;
    }
	
    /**
     * Returns an array of downlist-neighbor-atom lists.  Each list in the
     * array corresponds to a specific potential. A zero-length list indicates
     * that no concrete potentials apply to the atom.
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
    public void ensureCapacity(int index) {
        while (index > upList.length-1) {
            upList = (AtomArrayList[])Arrays.addObject(upList, new AtomArrayList());
            downList = (AtomArrayList[])Arrays.addObject(downList, new AtomArrayList());
        }
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
