/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLatticeNbrIterator;
import etomica.lattice.RectangularLatticeNbrIteratorAdjacent;
import etomica.potential.IteratorDirective;
import etomica.potential.IteratorDirective.Direction;

/**
 * Iteration is performed using site lists.
 * Direction is related to ordering of sites.
 */
public class Api1ASite implements AtomsetIteratorPDT, java.io.Serializable {
    
    /**
	 * Constructor makes iterator that must have box specified and then be 
	 * reset() before iteration.
     * 
     * @param D the dimension of the space of the simulation (used to construct cell iterators)
     */
	public Api1ASite(int D, BoxAgentManager agentManager) {
        neighborIterator = new RectangularLatticeNbrIteratorAdjacent(D);
        
        latticeIndex = new int[D];
        periodicity = new boolean[D];
        
        neighborIterator.setDirection(null);
        boxAgentManager = agentManager;
	}

	public void setBox(Box box) {
        neighborSiteManager = (NeighborSiteManager)boxAgentManager.getAgent(box);
        lattice = neighborSiteManager.getLattice();
        neighborIterator.setLattice(lattice);
        for (int i=0; i<periodicity.length; i++) {
            periodicity[i] = box.getBoundary().getPeriodicity(i);
        }
        neighborIterator.setPeriodicity(periodicity);
	}
    
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
	}
	
    public IAtomList next() {
        if (upListNow) {
            if (neighborIterator.hasNext()) {
                pair.atom1 = ((AtomSite)neighborIterator.next()).getAtom();
                return pair;
            }
            else if (!doGoDown) {
                return null;
            }
            upListNow = false;
            neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
            neighborIterator.reset();
            pair.atom1 = ((AtomSite)neighborIterator.next()).getAtom();
            return pair;
        }
        if (!neighborIterator.hasNext()) {
            return null;
        }
        pair.atom0 = ((AtomSite)neighborIterator.next()).getAtom();
        return pair;
    }

    public void unset() {
        neighborIterator.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    public void reset() {
        if(targetAtom == null) {
            unset();
            return;
        }

        lattice.latticeIndex(neighborSiteManager.getSite(targetAtom).latticeArrayIndex,latticeIndex);
        upListNow = (direction != IteratorDirective.Direction.DOWN);
        neighborIterator.setSite(latticeIndex);
        if (upListNow) {
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            pair.atom0 = targetAtom;
        }
        else {
            neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
            pair.atom1 = targetAtom;
        }
        neighborIterator.reset();
    }
    
    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. Direction is determined by ordering within occupant
     * list of cell of target atom, and then by the cell ordering of neighboring cells.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        doGoDown = (direction != IteratorDirective.Direction.UP);
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not 
     * in one of the species given at construction, no iterates will be returned.
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
    }

    
    private static final long serialVersionUID = 1L;
    private final RectangularLatticeNbrIterator neighborIterator;
    private final AtomPair pair = new AtomPair();
    protected final boolean[] periodicity;
    private final int[] latticeIndex;
    private boolean doGoDown;
    private boolean upListNow;
    private IteratorDirective.Direction direction;
    private IAtom targetAtom;
    private final BoxAgentManager boxAgentManager;
    private NeighborSiteManager neighborSiteManager;
    
    private CellLattice lattice;
    
}
