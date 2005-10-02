/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.list;

import java.io.IOException;
import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomTypeGroup;
import etomica.atom.SpeciesMaster;
import etomica.atom.SpeciesRoot;
import etomica.nbr.cell.AtomSequencerCell;
import etomica.util.Arrays;
import etomica.util.DirtyObject;
import etomica.util.EtomicaObjectInputStream;

/**
 * Sequencer used to maintain neighbor lists.  Holds lists of atoms
 * that were elsewhere deemed to be neighbors of the sequencer's atom.
 * Atoms are stored with reference to the potential that governs their
 * interactions.
 */
public class AtomSequencerNbr extends AtomSequencerCell implements DirtyObject {

    protected transient AtomArrayList[] upList, downList;
	
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
    public void ensureCapacity(int index) {
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
        System.out.println("writing");
        // write nothing.
        out.defaultWriteObject();
        
        int[][][] atomListInts = new int[2][][];
        atomListInts[0] = new int[upList.length][];
        for (int i=0; i<upList.length; i++) {
            atomListInts[0][i] = new int[upList[i].size()];
            for (int j=0; j<atomListInts[0].length; j++) {
                atomListInts[0][i][j] = upList[i].get(i).node.index();
                
            }
        }
        
        atomListInts[1] = new int[downList.length][];
        for (int i=0; i<upList.length; i++) {
            atomListInts[1][i] = new int[downList[i].size()];
            for (int j=0; j<atomListInts[1].length; j++) {
                atomListInts[1][i][j] = downList[i].get(i).node.index();
            }
        }
        ObjectData data = new ObjectData();
        data.atomListInts = atomListInts;
        // speciesAgent => SpeciesMaster
        data.speciesMaster = (SpeciesMaster)atom.node.parentSpeciesAgent().node.parentGroup();
        out.writeObject(data);
        System.out.println("done");
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
        ObjectData myData = (ObjectData)data;
        upList = new AtomArrayList[myData.atomListInts[0].length];
        for (int i=0; i<upList.length; i++) {
            for (int j=0; j<myData.atomListInts[0][i].length; j++) {
                upList[i].add(getAtomForIndex(myData.atomListInts[0][i][j],myData.speciesMaster));
            }
        }
        downList = new AtomArrayList[myData.atomListInts[1].length];
        for (int i=0; i<downList.length; i++) {
            for (int j=0; j<myData.atomListInts[1][i].length; j++) {
                upList[i].add(getAtomForIndex(myData.atomListInts[1][i][j],myData.speciesMaster));
            }
        }
    }
    
    /**
     * A factory class that will make this atom sequencer.  Typically this
     * is passed to the constructor of the atom factory.
     */
    public static final AtomSequencerFactory FACTORY = new AtomSequencerNbr.Factory();
    
    public static final class Factory implements AtomSequencerFactory, java.io.Serializable 
	{
        public AtomLinker makeSequencer(Atom atom) {
            return new AtomSequencerNbr(atom);
        }
    };

    /**
     * Wrapper class used to store both neighbor lists and SpeciesMaster
     */
    private static class ObjectData implements Serializable {
        int[][][] atomListInts;
        SpeciesMaster speciesMaster;
    }
    
    /**
     * Finds the child of the given speciesMaster having the given index.
     * @throws IllegalArgumentException if no Atom is found.
     */
    //Calling this repeatedly will be quite expensive for phases with long (flat)
    //AtomLists.  
    private static Atom getAtomForIndex(int index, SpeciesMaster speciesMaster) {
        Atom atom = speciesMaster;
outer:  while (true) {
            if (atom.node.index() == index) {
                return atom;
            }
            if (!(atom.type instanceof AtomTypeGroup)) {
                // atom of interest thinks it's below the leaf level
                throw new IllegalArgumentException("Failed to find Atom with index "+index);
            }
            AtomList childList = ((AtomTreeNodeGroup)atom.node).childList;
            for (AtomLinker link=childList.header.next; link!=childList.header; link=link.next) {
                if (link.atom.type.getIndexManager().sameAncestry(link.atom.node.index(),index)) {
                    atom = link.atom;
                    continue outer;
                }
            }
            
            // no children were ancestors of the atom of interest
            throw new IllegalArgumentException("Failed to find Atom with index "+index);
        }
    }
}
