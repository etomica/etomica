package etomica.atom.iterator;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Atom pair iterator that generates its iterates from a single basis atom, but
 * with inner- and outer-loop iterates taken from different sets of atoms
 * defined under the basis. Inner-loop iteration is "fixed" in the sense that
 * its iterates are independent of the outer-loop atom.
 * <p>
 * This iterator is constructed by extending ApiIntergroup, also forms atom
 * pairs by taking each atom from one of two disjoint sets of atoms. This feature
 * simplifies iteration becuase there needn't be concern for generating a pair
 * twice, or generating a "pair" that is just two handles to the same atom.
 * 
 * @author D.A Kofke and A.J. Schultz
 */
public class ApiIntragroupFixed extends ApiIntergroup implements
        AtomsetIteratorDirectable {

//do not make default iterator, because inner- and outer-loop iterates won't be disjoint sets 
    
    /**
     * Constructs a pair iterator that returns iterates from the given
     * pairIterator, which is expected to be a composite formed from
     * two basis-dependent iterators.
     */
    public ApiIntragroupFixed(ApiComposite pairIterator) {
        super(pairIterator);
    }

    /**
     * Indicates the direction of iteration to be performed with respect to
     * a specified target atom.  Has no effect if a target atom is not specified 
     * yet, or before iteration.
     */
    public void setDirection(Direction newDirection) {
        direction = newDirection;
        needSetupIterators = true;
    }

    protected void setupIterators() {
        if (direction == null || targetAtom == null) {
            // only need to worry about direction if there's a target
            aiOuter.setBasis(basisAtom);
            aiInner.setBasis(basisAtom);
            super.setupIterators();
        } else {
            if (aiInner.haveTarget(targetAtom)) {
                if (direction == IteratorDirective.Direction.DOWN) {
                    aiOuter.setBasis(basisAtom);
                    aiInner.setBasis(basisAtom);
                    aiOuter.setTarget(null);
                    aiInner.setTarget(targetAtom);
                } else {
                    // can't iterate "up" if the target type matches the
                    // first iterator.
                    aiOuter.setBasis(null);
                }
            } else {
                if (direction == IteratorDirective.Direction.UP) {
                    aiOuter.setBasis(basisAtom);
                    aiInner.setBasis(basisAtom);
                    aiOuter.setTarget(targetAtom);
                    aiInner.setTarget(null);
                } else {
                    // can't iterate "down" if the target type matches the
                    // second iterator
                    aiOuter.setBasis(null);
                }
            }
            needSetupIterators = false;
        }
    }

    /**
     * Specifies the basis, which identifies the atoms subject to iteration. The
     * given atomSet must be of length 1, otherwise the iterator is unset. The
     * atom specifies the basis for both the outer-loop and inner-loop
     * iteration. In each case, if the basis atom is not a leaf atom, its
     * children will be the subject of iteration. If the basis atom is a leaf,
     * it will itself be the iterate. If given atomset is null, or if its length
     * is not equal to 1, iterator will give no iterates until a proper basis is
     * specified via another call to this method.
     */
    public void setBasis(AtomSet basisAtoms) {
        if (basisAtoms == null || basisAtoms.getAtomCount() != 1) {
            basisAtom = null;
        } else {
            basisAtom = basisAtoms.getAtom(0);
        }
        needSetupIterators = true;
    }

    /**
     * Returns 1, indicating that the basis should specify one atom.
     */
    public int basisSize() {
        return 1;
    }

    protected IAtom basisAtom;
    protected Direction direction;
    private static final long serialVersionUID = 1L;
}
