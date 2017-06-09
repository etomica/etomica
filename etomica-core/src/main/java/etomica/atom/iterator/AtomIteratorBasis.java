/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.AtomArrayList;

/**
 * Elementary basis-dependent iterator that gives atoms meeting specification
 * of a basis and a target.  Iterates are determined as follows:
 * <ul>
 * 
 * <li>the <i>basis atom</i> determines the set of atoms (basis-set) that are candidates for iteration.
 * If the basis is an atom group, iterates (if any) will always be the child 
 * atoms of the basis; if the basis is a leaf atom, the single iterate (if any)
 * will be the basis atom itself.
 * 
 * <li>the <i>target atom</i> narrows down the iterates from the basis-set. Possibilities
 * are: 
 * <ul>
 * <li>no target atom is specified: all basis-set atoms will be given on iteration
 * <li>the target atom is the same as the basis atom: this is equivalent to specifying no target
 * <li>the target atom is in the species hierarchy above the basis atom, and the
 * basis atom is descended from it: this is equivalent to specifying no target
 * <li>the target atom is a child of the basis atom: the target atom will be the only iterate
 * <li>the target atom is descended from the basis atom, but is not a direct child of it:
 * the only iterate will be the child atom of the basis through which the target is descended
 * <li>the target and basis are in different branches of the species hierarchy, such that
 * the target atom is not descended from the basis atom, and the basis atom is not descended
 * from the target atom: no iterates will be given
 * </ul>
 *
 * </ul>
 */
public class AtomIteratorBasis extends AtomIteratorArrayListSimple implements
        AtomIteratorBasisDependent {

    /**
     * Constructor makes iterator in an unset condition; must set basis and call reset before
     * beginning iteration.
     */
    public AtomIteratorBasis() {
        super();
        littleList.clear();
        myList = littleList;
    }
    
    

    /**
     * Method to specify a target atom. Specifying a zero-length AtomSet or a
     * length-1 AtomSet with a null atom releases
     * any target restrictions, and specifies that the iterator should give all
     * of the basis-set atoms. Call to this method leaves iterator unset; call to reset is
     * required before beginning iteration.
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
        needSetupIterator = (basis != null);//flag to setup iterator only if
                                            // presently has a non-null basis
        unset();
    }

    /**
     * Sets the basis for iteration, such that the childList atoms of the given
     * atom will be subject to iteration (within any specifications given by a
     * prior or subsequent call to setTarget). If given atom is a leaf, it will
     * itself be the sole candidate iterate given by the iterator. If argument is null or
     * otherwise does not specify an atom, iterator will be conditioned to give
     * no iterates until a new basis is specified. The given AtomSet, if not
     * null, must have a size of 1.
     * 
     * @throws IllegalArgumentException 
     *              if atoms.count() is not 0 or 1
     */
    public void setBasis(IMoleculeList atoms) {
        if (atoms == null) {
            basis = null;
            littleList.clear();
            myList = littleList;
            setList(myList);
            needSetupIterator = false;
        } else if (atoms.getMoleculeCount() == 1) {
            basis = atoms.getMolecule(0);
            needSetupIterator = true;
        } else {
            throw new IllegalArgumentException(
                    "Inappropriate number of atoms given in basis");
        }
        unset();
    }

    /**
     * Returns true if the given target with the present basis could
     * yield an iterate. Assumes that the basis -- if it is a group -- 
     * has child atoms. 
     */
    public boolean haveTarget(IAtom target) {
        if(basis == null) return false;
        if (target == null) {
            return true;
        }
        return target.getParentGroup() == basis;
    }

    /**
     * Puts iterator in a state ready to begin iteration.
     */
    public void reset() {
        if (basis == null) {
            return;
        }
        if (needSetupIterator) {
            setupIterator();
        }
        setList(myList);
        super.reset();
    }

    /**
     * Returns 1, indicating that only a single-atom basis is appropriate.
     */
    public int basisSize() {
        return 1;
    }

    /**
     * Common method to complete tasks needed to adjust to new target or basis.
     * Any call to setBasis or setTarget sets flag that indicates this method
     * should be invoked upon reset.
     */
    private void setupIterator() {
        needSetupIterator = false;
        littleList.clear();
        myList = littleList;
        if (basis != null) {
            if (targetAtom == null) {
                setupBasisIteration();
            }
            else if (targetAtom.getParentGroup() == basis) {
                //targetAtom is the child of the basis atom
                littleList.add(targetAtom);
            }
        }
    }

    /**
     * Convenience method used by setupIterator
     */
    private void setupBasisIteration() {
        myList = basis.getChildList();
    }

    private static final long serialVersionUID = 1L;
    private final AtomArrayList littleList = new AtomArrayList(1);//used to form a list of
                                                       // one iterate if target
                                                       // is specified
    private IAtom targetAtom;
    private IMolecule basis;
    private IAtomList myList;
    private boolean needSetupIterator = true;//flag to indicate if
                                             // setupIterator must be called
                                             // upon reset
}
