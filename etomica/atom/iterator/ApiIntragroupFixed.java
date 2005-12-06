package etomica.atom.iterator;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.iterator.IteratorDirective.Direction;

public class ApiIntragroupFixed extends ApiIntergroup implements
        AtomsetIteratorDirectable {

    public ApiIntragroupFixed() {
        super();
    }

    public ApiIntragroupFixed(ApiComposite pairIterator) {
        super(pairIterator);
    }

    public void setDirection(Direction newDirection) {
        direction = newDirection;
        needSetupIterators = true;
    }
    
    protected void setupIterators() {
        if (direction == null || targetAtoms.count() != 1) {
            // only need to worry about direction if there's 1 target
            aiOuter.setBasis(basisAtom);
            aiInner.setBasis(basisAtom);
            super.setupIterators();
        }
        else {
            Atom target = targetAtoms.getAtom(0);
            if (aiInner.haveTarget(target)) {
                if (direction == IteratorDirective.DOWN) {
                    aiOuter.setBasis(basisAtom);
                    aiInner.setBasis(basisAtom);
                    aiOuter.setTarget(emptyTarget);
                    aiInner.setTarget(target);
                }
                else {
                    // can't iterate "down" if the target type matches the
                    // first iterator.
                    aiOuter.setBasis(null);
                }
            } else {
                if (direction == IteratorDirective.UP) {
                    aiOuter.setBasis(basisAtom);
                    aiInner.setBasis(basisAtom);
                    aiOuter.setTarget(target);
                    aiInner.setTarget(emptyTarget);
                }
                else {
                    // can't iterate "up" if the target type matches the
                    // second iterator
                    aiOuter.setBasis(null);
                }
            }
            needSetupIterators = false;
        }
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
        if (basisAtoms == null || basisAtoms.count() != 1) {
            basisAtom = null;
        } else {
            basisAtom = basisAtoms.getAtom(0);
        }
        needSetupIterators = true;
    }

    public int basisSize() {
        return 1;
    }
    
    protected Atom basisAtom;
    protected Direction direction;
}
