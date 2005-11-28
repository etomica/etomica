/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomsetArray;

/**
 * Iterator that returns pairs formed using two different basis atoms, so that
 * the iterates are taken from two different groups.
 */
public final class ApiIntergroup extends AtomPairIteratorAdapter implements
        AtomsetIteratorBasisDependent, ApiComposite {

    public ApiIntergroup() {
        this(new ApiInnerFixed(new AtomIteratorBasis(),
                  new AtomIteratorBasis()));
    }

    /**
     * Constructs a pair iterator that returns iterates from the given
     * pairIterator, which is expected to contain two basis-dependent 
     * iterators.
     */
    public ApiIntergroup(ApiComposite pairIterator) {
        super(pairIterator);
        aiOuter = (AtomsetIteratorBasisDependent) pairIterator
                .getOuterIterator();
        aiInner = (AtomsetIteratorBasisDependent) pairIterator
                .getInnerIterator();
    }

    /**
     * @throws NullPointerException
     *          if targetAtoms is null; use AtomSet.NULL instead
     * @throws IllegalArgumentException
     *          if targetAtoms.count() > 2
     */
    public void setTarget(AtomSet targetAtoms) {
        if (targetAtoms == null)
            throw new NullPointerException (
                    "Cannot set target to null; use AtomSet.NULL");
        if (targetAtoms.count() > 2) {
            throw new IllegalArgumentException(
                    "Too many target atoms for iterator");
        }
        this.targetAtoms = targetAtoms;
        needSetupIterators = true;
    }

    public boolean haveTarget(AtomSet targets) {
        switch (targets.count()) {
        case 0:
            return true;
        case 1:
            Atom target = targets.getAtom(0);
            return aiOuter.haveTarget(target) || aiInner.haveTarget(target);
        case 2:
            Atom target0 = targets.getAtom(0);
            Atom target1 = targets.getAtom(1);
            return (aiOuter.haveTarget(target0) && aiInner.haveTarget(target1))
                    || (aiOuter.haveTarget(target1) && aiInner
                            .haveTarget(target0));
        default:
            return false;
        }
    }

    private void setupIterators() {
        switch (targetAtoms.count()) {
        case 0:
            aiOuter.setTarget(targetAtoms);
            aiInner.setTarget(targetAtoms);
            break;
        case 1:
            Atom target = targetAtoms.getAtom(0);
            if (aiInner.haveTarget(target)) {
                aiOuter.setTarget(emptyTarget);
                aiInner.setTarget(target);
            } else {
                aiOuter.setTarget(target);
                aiInner.setTarget(emptyTarget);
            }
            break;
        case 2:
            Atom target0 = targetAtoms.getAtom(0);
            Atom target1 = targetAtoms.getAtom(1);
            if (aiInner.haveTarget(target0)) {
                aiOuter.setTarget(target1);
                aiInner.setTarget(target0);
            } else {
                aiOuter.setTarget(target0);
                aiInner.setTarget(target1);
            }
            break;
        default:
            throw new IllegalArgumentException(
                    "Too many target atoms for iterator");
        }
        needSetupIterators = false;
    }

    /**
     * Specifies the basis, which identifies the atoms subject to iteration. The
     * given atomSet must be of length 2. The first atom in the set specifies
     * the basis for the outer-loop iteration, and second atom specifies the
     * basis for the inner-loop iteration. In each case, if the basis atom is
     * not a leaf atom, its children will be the subject of iteration. If the
     * basis atom is a leaf, it will itself be the iterate. If given atomset is
     * null, or if its length is not equal to 2, iterator will give no iterates
     * until a proper basis is specified via another call to this method.
     */
    public void setBasis(AtomSet basisAtoms) {
        if (basisAtoms == null || basisAtoms.count() != 2) {
            aiOuter.setBasis(null);
        } else {
            aiOuter.setBasis(basisAtoms.getAtom(0));
            aiInner.setBasis(basisAtoms.getAtom(1));
        }
        needSetupIterators = true;
    }

    public void reset() {
        if (needSetupIterators) {
            setupIterators();
        }
        super.reset();
    }

    /**
     * Returns 2, indicating that the setBasis method expects an array of two
     * atoms.
     */
    public final int basisSize() {
        return 2;
    }

    /**
     * Returns the (basis-dependent) iterator instance used for the inner-loop
     * iteration
     */
    public AtomIterator getInnerIterator() {
        return (AtomIterator) aiInner;
    }

    /**
     * Returns the (basis-dependent) iterator instance used for the outer-loop
     * iteration
     */
    public AtomIterator getOuterIterator() {
        return (AtomIterator) aiOuter;
    }
    
    /**
     * Sets up atom iterators (if needed) and performs action on all 
     * pair iterates.
     */
    public void allAtoms(AtomsetAction action) {
        if(needSetupIterators) {
            setupIterators();
        }
        super.allAtoms(action);
    }

    private final AtomsetIteratorBasisDependent aiOuter;
    private final AtomsetIteratorBasisDependent aiInner;
    private AtomSet targetAtoms = AtomSet.NULL;
    private final AtomsetArray emptyTarget = new AtomsetArray(0);
    private boolean needSetupIterators = true;

}