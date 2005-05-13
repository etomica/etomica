package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.atom.AtomTreeNode;
import etomica.atom.AtomsetArray;


/**
 * Singlet pair iterator for which atoms in pair are defined via
 * a basis pair and a target specification.  Name indicates that it returns
 * 1 atom from one basis, and 1 atom from the other.  Basis may be same or different
 * from each other.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 10, 2005 by kofke
 */
public class Api11 extends AtomPairIteratorAdapter implements
        AtomsetIteratorBasisDependent {

    /**
     * Creates iterator in state not ready for iteration.
     * Must set basis, target and then reset before iterating.
     */
    public Api11() {
        super(new ApiSinglet());
        needUpdateAtoms = true;
        canReset = false;
    }

    /**
     * Defines the basis, which must have two different atoms.  The atoms in the pair
     * iterate will be those child atoms through which the specified target
     * atoms are descended.
     */
    public void setBasis(AtomSet atoms) {
        if(atoms.count() != 2) throw new IllegalArgumentException("basis must have two atoms");
        basisNode0 = atoms.getAtom(0).node;
        basisNode1 = atoms.getAtom(1).node;
        needUpdateAtoms = true;
    }

    /**
     * Returns 2, indicating that the basis takes two atoms.
     */
    public int basisSize() {
        return 2;
    }

    /**
     * Returns true if the indicated target atoms would yield an iterate
     * for the current basis.
     */
    public boolean haveTarget(AtomSet newTarget) {
        if(newTarget.count() < 2) return false;
        //save current target
        if(oldTarget.count() != targetAtoms.count()) {
            oldTarget = new AtomsetArray(targetAtoms);
        } else {
            oldTarget.setAtoms(targetAtoms);
        }
        boolean hasNext = iterator.hasNext();
        //establish if haveTarget by going through setTarget process
        setTarget(newTarget);
        updateAtoms();
        boolean haveTarget = canReset;
        //then restore old target
        targetAtoms.setAtoms(oldTarget);
        needUpdateAtoms = true;
        if(hasNext) reset();
        return haveTarget;
    }

    /**
     * Sets target atoms for iteration. Atom set must contain at least two
     * atoms, and there must be at least one atom descended from each basis
     * atom.  Must identify exactly two child atoms (one for each basis) through
     * which all target atoms are descended.  If these conditions are not met,
     * then no iterate will be given on reset.
     */
    public void setTarget(AtomSet newTarget) {
        if(newTarget.count() < 2) throw new IllegalArgumentException("11 iterator requires at least two target atoms");
        if(targetAtoms.count() != newTarget.count()) {
            targetAtoms = new AtomsetArray(newTarget);
        } else {
            targetAtoms.setAtoms(newTarget);
        }   
        needUpdateAtoms = true;
    }
    
    private void updateAtoms() {
        if(basisNode0 == null || basisNode1 == null) {
            canReset = false;
            return;
        }
        AtomTreeNode target0 = null;
        AtomTreeNode target1 = null;
        canReset = true;
        for(int i=targetAtoms.count()-1; i>=0; i--) {
            Atom target = targetAtoms.getAtom(i);
            if(target == null) continue;
            //see if target specifies basis0 child
            AtomTreeNode newTarget = target.node.childWhereDescendedFrom(basisNode0);
            if(newTarget != null) {//it does
                if(target0 == null) {//found basis0 child
                    target0 = newTarget;
                } else if(target0 != newTarget) {//target specifies multiple atoms under basis0
                    if(basisNode0 != basisNode1) {//that's bad if basis1 is not the same as basis0
                        canReset = false;
                        break;
                    } else if(target1 == null) {//otherwise assign target1 if not assigned before
                        target1 = newTarget;
                    } else if(target1 != newTarget) {
                        canReset = false;
                        break;
                    }
                }
            } else {//it doesn't; try for basis1
                newTarget = target.node.childWhereDescendedFrom(basisNode1);
                if(newTarget != null) {//it does
                    if(target1 == null) {//found basis1 child
                        target1 = newTarget;
                    } else if(target1 != newTarget) {//target give multiple atoms under basis1
                        canReset = false;
                        break;
                    }
                } else {//target atom not derived from either basis
                    canReset = false;
                    break;
                }
            }
        }//end for
        if(target0 == null || target1 == null) {
            canReset = false;
        } else {
            ((ApiSinglet)iterator).setPair(target0.atom(), target1.atom());
        }
        needUpdateAtoms = false;
    }
    
    /**
     * Resets for iteration if an atom pair has been adequately 
     * specified via previous calls to setBasis and setTarget.
     * Otherwise iterator is unset.
     */
    public void reset() {
        if(needUpdateAtoms) updateAtoms();
        if(canReset) {
            super.reset();
        } else {
            unset();
        }
    }
    
    /**
     * Performs action on atom pair if one has been adequately
     * specified via previous calls to setBasis and setTarget.
     */
    public void allAtoms(AtomsetAction action) {
        if(canReset) super.allAtoms(action);
    }
    
    private AtomTreeNode basisNode0, basisNode1;
    private AtomsetArray targetAtoms = new AtomsetArray(0);
    private AtomsetArray oldTarget = new AtomsetArray(0);
    private boolean needUpdateAtoms;
    private boolean canReset;
}
