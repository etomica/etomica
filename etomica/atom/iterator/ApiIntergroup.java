package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;

/**
 * Iterator that returns pairs formed using two different basis atoms, so that
 * the iterates are taken from two different groups.
 */
public class ApiIntergroup extends AtomPairIteratorAdapter implements
        AtomsetIteratorBasisDependent, ApiComposite {

    /**
     * Default iterator is an ApiInnerFixed formed from two AtomIteratorBasis
     * instances.
     */
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
     * Specifies a target atom, which should appear in all iterates. A
     * null value removes any restriction on the iterates.
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
        needSetupIterators = true;
    }

    public boolean haveTarget(IAtom target) {
        if (target == null) {
            return true;
        }
        return aiOuter.haveTarget(target) || aiInner.haveTarget(target);
    }

    protected void setupIterators() {
        if (targetAtom == null) {
            aiOuter.setTarget(targetAtom);
            aiInner.setTarget(targetAtom);
        }
        else {
            if (aiInner.haveTarget(targetAtom)) {
                aiOuter.setTarget(null);
                aiInner.setTarget(targetAtom);
            } else {
                aiOuter.setTarget(targetAtom);
                aiInner.setTarget(null);
            }
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
        if (basisAtoms == null || basisAtoms.getAtomCount() != 2) {
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
    public int basisSize() {
        return 2;
    }

    /**
     * Returns the (basis-dependent) iterator instance used for the inner-loop
     * iteration.
     */
    public AtomIterator getInnerIterator() {
        return (AtomIterator) aiInner;
    }

    /**
     * Returns the (basis-dependent) iterator instance used for the outer-loop
     * iteration.
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

    private static final long serialVersionUID = 1L;
    protected final AtomsetIteratorBasisDependent aiOuter;
    protected final AtomsetIteratorBasisDependent aiInner;
    protected IAtom targetAtom;
    protected boolean needSetupIterators = true;

}