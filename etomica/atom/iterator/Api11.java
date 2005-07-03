package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.AtomTreeNode;
import etomica.action.AtomsetAction;
import etomica.atom.AtomsetArray;

/**
 * Singlet pair iterator for which atoms in pair are defined via a basis pair
 * and a target specification. Name of class indicates that it returns 1 atom
 * from one basis, and 1 atom from the other. Basis may be same or different
 * from each other. Atoms specified by basis and target are largely as defined
 * in AtomIteratorBasis, but with some differences. Unlike AtomIteratorBasis, no
 * iterate atom is indicated if target doesn't specify it, so 
 * <ul>
 * <li>target atom
 * set must include at least two atoms (anything less throws an exception); 
 * <li>more than two target atoms can be specified, but they must consistently point
 * to exactly two iterate atoms; and 
 * <li>presence of any target atom that
 * doesn't indicate an atom for the current basis causes iterator to yield no
 * iterates (this contrasts with AtomIteratorBasis, for which a target atom
 * could be in hierarchy above the basis atom). 
 * </ul>
 * Expected use of this class is inside another iterator that directs targets of
 * fewer than two atoms to 1A or AA iterators. <p>
 * Rules for determining atoms in iterate pair are as follows.
 * <ul>
 * <li>The basis locates the candidate atoms, with each atom in the iterate
 * pair determined by one of the basis atoms. The first atom in the iterate pair
 * corresponds to the first basis atom, and the second atom in the iterate pair
 * corresponds to the second basis atom. In each case, the candidate iterate
 * atom is a child of the basis atom, unless the basis is a leaf atom, in which
 * case the candidate iterate atom will be the basis atom itself. If basis atoms
 * are the same as each other, the iterate atoms will be ordered such that
 * atom0.compareTo(atom1) < 0, otherwise they are ordered according to the basis
 * atom ordering.
 * <li>The target specifies which basis child-atom appears in the iterate pair.
 * There must be at least two non-null target atoms given in the target atom
 * set. The order of the atoms in the target atom set is irrelevant. The target
 * must specify exactly one atom for each basis atom (or two if both basis atoms
 * are the same atom). If either of the basis atoms has no atoms corresponding
 * to it, or has multiple atoms (if basis atoms are different) corresponding to
 * it, then no iterates are given.
 * <li>Every target atom must specify an iterate atom, so if there are more
 * than two target atoms, some target atoms must specify the same iterate atom
 * (e.g., they may be different descendants of the same child of a basis atom).
 * <li>A target atom will indicate an iterate atom for a basis if the target is
 * descended from the basis; then the iterate atom will be the child of the
 * basis through which the target atom is descended (which is the target itself
 * if it is a child of the basis atom). If the basis is a leaf atom, then the
 * target specifies an atom only if it is the basis atom itself, and then the
 * iterate will be the basis/target atom.
 * </ul>
 * 
 * 
 * @author David Kofke
 * @see AtomIteratorBasis
 *  
 */

/*
 * History Created on May 10, 2005 by kofke
 */
public class Api11 extends AtomPairIteratorAdapter implements
        AtomsetIteratorBasisDependent {

    /**
     * Creates iterator in state not ready for iteration. Must set basis, target
     * and then reset before iterating.
     */
    public Api11() {
        super(new ApiSinglet());
        needUpdateAtoms = true;
        canReset = false;
    }

    /**
     * Defines the basis, which must have two different atoms. The atoms in the
     * pair iterate will be those child atoms through which the specified target
     * atoms are descended, or the basis atom itself if the basis is a leaf atom
     * and consistent with the target.
     */
    public void setBasis(AtomSet atoms) {
        if (atoms.count() != 2)
            throw new IllegalArgumentException("basis must have two atoms");
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
     * Returns true if the indicated target atoms would yield an iterate for the
     * current basis. Returns false if no iterate would be given or if
     * newTarget.count() < 2.
     */
    public boolean haveTarget(AtomSet newTarget) {
        if (newTarget.count() < 2) {
            return false;
        }
        //save current target
        if (oldTarget.count() != targetAtoms.count()) {
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
        if (oldTarget.count() != targetAtoms.count()) {
            targetAtoms = new AtomsetArray(oldTarget);
        } else {
            targetAtoms.setAtoms(oldTarget);
        }
        needUpdateAtoms = true;
        if (hasNext) {
            reset();
        }
        return haveTarget;
    }

    /**
     * Sets target atoms for iteration. Atom set must contain at least two
     * atoms. Target must identify exactly two child atoms (one for each basis)
     * through which all target atoms are descended. If these conditions are not
     * met, then no iterate will be given on reset.
     * 
     * @throws IllegalArgumentException
     *             if newTarget.count() < 2
     * @throws NullPointerException 
     *             if newTarget is null
     */
    public void setTarget(AtomSet newTarget) {
        if (newTarget.count() < 2)
            throw new IllegalArgumentException(
                    "11 iterator requires at least two target atoms");
        if (targetAtoms.count() != newTarget.count()) {
            targetAtoms = new AtomsetArray(newTarget);
        } else {
            targetAtoms.setAtoms(newTarget);
        }
        needUpdateAtoms = true;
    }

    private void updateAtoms() {
        if (basisNode0 == null || basisNode1 == null) {
            canReset = false;
            return;
        }
        AtomTreeNode target0 = null;
        AtomTreeNode target1 = null;
        canReset = true;
        for (int i = targetAtoms.count() - 1; i >= 0; i--) {
            Atom target = targetAtoms.getAtom(i);
            if (target == null)
                continue;
            //see if target specifies basis0 child
            AtomTreeNode newTarget = null;
            if (basisNode0.isLeaf() && target.node == basisNode0)
                newTarget = basisNode0;
            else
                newTarget = target.node.childWhereDescendedFrom(basisNode0);
            if (newTarget != null) {//it does
                if (target0 == null) {//found basis0 child
                    target0 = newTarget;
                } else if (target0 != newTarget) {//target specifies multiple
                    // atoms under basis0
                    if (basisNode0 != basisNode1) {//that's bad if basis1 is
                        // not the same as basis0
                        canReset = false;
                        break;
                    } else if (target1 == null) {//otherwise assign target1 if
                        // not assigned before
                        if (newTarget.compareTo(target0) > 0) {
                            target1 = newTarget;
                        } else {
                            target1 = target0;
                            target0 = newTarget;
                        }
                    } else if (target1 != newTarget) {
                        canReset = false;
                        break;
                    }
                }
            } else {//it doesn't; try for basis1
                if (basisNode1.isLeaf() && target.node == basisNode1)
                    newTarget = basisNode1;
                else
                    newTarget = target.node.childWhereDescendedFrom(basisNode1);
                if (newTarget != null) {//it does
                    if (target1 == null) {//found basis1 child
                        target1 = newTarget;
                    } else if (target1 != newTarget) {//target give multiple
                        // atoms under basis1
                        canReset = false;
                        break;
                    }
                } else {//target atom not derived from either basis
                    canReset = false;
                    break;
                }
            }
        }//end for
        if (target0 == null || target1 == null) {
            canReset = false;
        } else {
            ((ApiSinglet) iterator).setPair(target0.atom(), target1.atom());
        }
        needUpdateAtoms = false;
    }

    /**
     * Resets for iteration if an atom pair has been adequately specified via
     * previous calls to setBasis and setTarget. Otherwise iterator is unset.
     */
    public void reset() {
        if (needUpdateAtoms)
            updateAtoms();
        if (canReset) {
            super.reset();
        } else {
            unset();
        }
    }

    /**
     * Returns 1 if the current basis and target will give an iterate, otherwise
     * returns 0.
     */
    public int size() {
        if (needUpdateAtoms)
            updateAtoms();
        return canReset ? 1 : 0;
    }

    /**
     * Returns true if the current basis and target would give the indicate set
     * of atoms.
     */
    public boolean contains(AtomSet atoms) {
        if (needUpdateAtoms)
            updateAtoms();
        return canReset ? super.contains(atoms) : false;
    }

    /**
     * Performs action on atom pair if one has been adequately specified via
     * previous calls to setBasis and setTarget.
     */
    public void allAtoms(AtomsetAction action) {
        if (needUpdateAtoms)
            updateAtoms();
        if (canReset)
            super.allAtoms(action);
    }

    private AtomTreeNode basisNode0, basisNode1;
    private AtomsetArray targetAtoms = new AtomsetArray(0);
    private AtomsetArray oldTarget = new AtomsetArray(0);
    private boolean needUpdateAtoms;
    private boolean canReset;
}