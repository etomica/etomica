/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;

/**
 * Iterator that returns pairs formed using two different basis atoms, so that
 * the iterates are taken from two different groups.
 */
public final class ApiIntergroup extends AtomsetIteratorAdapter implements
        AtomsetIteratorBasisDependent, ApiComposite {

    public ApiIntergroup() {
        this(
                new ApiInnerFixed(new AtomIteratorBasis(),
                        new AtomIteratorBasis()));
    }

    public ApiIntergroup(ApiComposite pairIterator) {
        super(pairIterator);
        aiOuter = (AtomsetIteratorBasisDependent) pairIterator
                .getOuterIterator();
        aiInner = (AtomsetIteratorBasisDependent) pairIterator
                .getInnerIterator();
    }

    /*
     * (non-Javadoc)
     * 
     * @see etomica.AtomsetIteratorBasisDependent#setDirective(etomica.IteratorDirective)
     */
    public void setTarget(Atom[] targetAtoms) {
        aiOuter.setTarget(targetAtoms);
    }

    /**
     * Specifies the basis, which identifies the atoms subject to iteration. The
     * given array should be of length 2 (at least); first atom in array
     * specifies the basis for the outer-loop iteration, and second atom
     * specifies the basis for the inner-loop iteration. In each case, if the
     * basis atom is not a leaf atom, its children will be the subject of
     * iteration. If the basis atom is a leaf, it will itself be the iterate. If
     * given array is null, or if its length is not at least 2, iterator will
     * give no iterates until a proper basis is specified via another call to
     * this method.
     */
    public void setBasis(Atom[] basisAtoms) {
        if (basisAtoms == null || basisAtoms.length < 2) {
            aiOuter.setBasis(null);
        } else {
            atom[0] = basisAtoms[0];
            aiOuter.setBasis(atom);
            atom[0] = basisAtoms[1];
            aiInner.setBasis(atom);
        }
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

    private final AtomsetIteratorBasisDependent aiOuter;
    private final AtomsetIteratorBasisDependent aiInner;
    private final Atom[] atom = new Atom[1];

}