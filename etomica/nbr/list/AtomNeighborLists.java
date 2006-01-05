/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.list;

import java.io.IOException;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.util.Arrays;
import etomica.util.DirtyObject;
import etomica.util.EtomicaObjectInputStream;

/**
 * Sequencer used to maintain neighbor lists.  Holds lists of atoms
 * that were elsewhere deemed to be neighbors of the sequencer's atom.
 * Atoms are stored with reference to the potential that governs their
 * interactions.
 */
public class AtomNeighborLists implements DirtyObject {

    protected transient AtomArrayList[] upList, downList;
	
    /**
     * Constructs sequencer for the given atom.
     */
    public AtomNeighborLists() {
        upList = new AtomArrayList[0];
        downList = new AtomArrayList[0];
    }
    
    /**
     * Adds the given atom to the neighbor set that are uplist of
     * this sequencer's atom and which interact via the given potential.  
     * @param a the new uplist neighbor atom
     * @param potential the potential between the atoms
     */
    public void addUpNbr(Atom a, int index) {
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
    public void addDownNbr(Atom a, int index) {
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
    protected void ensureCapacity(int index) {
        while (index > upList.length-1) {
            upList = (AtomArrayList[])Arrays.addObject(upList, new AtomArrayList());
            downList = (AtomArrayList[])Arrays.addObject(downList, new AtomArrayList());
        }
    }
	
    /**
     * Should be called when removing a potential that applied to this atom
     */
    //TODO consider whether this method might foul up the index assignments
	public void decrementNbrListArrays() {
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
     * Write out neighbor lists as arrays of Atom indices.  Include
     * the SpeciesMaster so that rebuild() can find the Atoms corresponding
     * to those indices.
     */
    private void writeObject(java.io.ObjectOutputStream out)
    throws IOException
    {
        // write nothing.
        out.defaultWriteObject();
        
        int[][][] atomListInts = new int[2][][];
        atomListInts[0] = new int[upList.length][];
        for (int i=0; i<upList.length; i++) {
            atomListInts[0][i] = new int[upList[i].size()];
            for (int j=0; j<atomListInts[0][i].length; j++) {
                atomListInts[0][i][j] = upList[i].get(j).node.index();
                
            }
        }
        
        atomListInts[1] = new int[downList.length][];
        for (int i=0; i<upList.length; i++) {
            atomListInts[1][i] = new int[downList[i].size()];
            for (int j=0; j<atomListInts[1][i].length; j++) {
                atomListInts[1][i][j] = downList[i].get(j).node.index();
            }
        }
        out.writeObject(atomListInts);
    }

    /**
     * Just reads the data from the stream and stashes it back into the stream
     * object for later use.
     */
    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        EtomicaObjectInputStream etomicaIn = (EtomicaObjectInputStream)in; 
        //read nothing.
        etomicaIn.defaultReadObject();

        etomicaIn.objectData.put(this,etomicaIn.readObject());
        etomicaIn.dirtyObjects.add(this);
    }

    /**
     * Rebuilds the neighbor lists from the previously-read Atom indices.
     */
    public void rebuild(Object data) {
        int[][][] atomListInts = (int[][][])data;
        upList = new AtomArrayList[atomListInts[0].length];
        for (int i=0; i<upList.length; i++) {
            upList[i] = new AtomArrayList();
            for (int j=0; j<atomListInts[0][i].length; j++) {
                upList[i].add(EtomicaObjectInputStream.getAtomForIndex(atomListInts[0][i][j]));
            }
        }
        downList = new AtomArrayList[atomListInts[1].length];
        for (int i=0; i<downList.length; i++) {
            downList[i] = new AtomArrayList();
            for (int j=0; j<atomListInts[1][i].length; j++) {
                upList[i].add(EtomicaObjectInputStream.getAtomForIndex(atomListInts[1][i][j]));
            }
        }
    }
    
}
