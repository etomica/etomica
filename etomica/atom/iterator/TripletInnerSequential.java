package etomica.atom.iterator;
import java.io.Serializable;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.AtomSet;
import etomica.atom.AtomsetArray;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.IteratorDirective.Direction;

public class TripletInnerSequential implements AtomsetIteratorBasisDependent, 
        AtomsetIteratorDirectable, Serializable  {

    public TripletInnerSequential() {
        super();
        next = new AtomsetArray(3);
        atomArray = next.getArray();
    }

    public void setDirection(Direction direction) {
        doGoDown = true;
        doGoUp = true;
        if (direction == IteratorDirective.Direction.UP) {
            doGoDown = false;
        }
        else if (direction == IteratorDirective.Direction.DOWN) {
            doGoUp = false;
        }
    }

    public void setBasis(AtomSet atoms) {
        childList = ((IAtomGroup)atoms).getChildList();
    }

    public int basisSize() {
        return 1;
    }

    public boolean haveTarget(IAtom target) {
        int n = childList.getAtomCount();
        for (int i=0; i<n; i++) {
            if (childList.getAtom(i) == targetAtom) {
                return true;
            }
        }
        return false;
    }

    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
    }

    public void reset() {
        if (targetAtom != null) {
            int n = childList.getAtomCount();
            cursor = -1;
            for (int i=0; i<n; i++) {
                if (childList.getAtom(i) == targetAtom) {
                    cursor = i;
                    break;
                }
            }
            if (cursor == -1) {
                stateUpDown = 3;
            }
            else if (doGoUp) {
                if (cursor < childList.getAtomCount() - 2) {
                    stateUpDown = 0;
                }
                else if (doGoDown) {
                    if (cursor > 0 && cursor < childList.getAtomCount() - 1) {
                        stateUpDown = 1;
                    }
                    else if (cursor > 1 && cursor < childList.getAtomCount()) {
                        stateUpDown = 2;
                    }
                    else {
                        stateUpDown = 3;
                    }
                }
            }
            else {
                if (cursor > 1 && cursor < childList.getAtomCount()) {
                    stateUpDown = 2;
                }
                else {
                    stateUpDown = 3;
                }
            }
        }
        else {
            doGoUp = true;
            doGoDown = false;
            stateUpDown = 0;
            cursor = 0;
        }
    }

    public void unset() {
        cursor = childList.getAtomCount();
    }

    public AtomSet next() {
        atomArray[0] = childList.getAtom(cursor-stateUpDown);
        atomArray[1] = childList.getAtom(cursor-stateUpDown+1);
        atomArray[2] = childList.getAtom(cursor-stateUpDown+2);
        if (stateUpDown == 0) {
            if (targetAtom == null) {
                cursor++;
            }
            else if (doGoDown) {
                stateUpDown = 1;
            }
        }
        else if (stateUpDown == 1) {
            stateUpDown = 2;
        }
        else if (stateUpDown == 2) {
            stateUpDown = 3;
        }
        return next;
    }

    public void allAtoms(AtomsetAction action) {
        reset();
        for (AtomSet atoms = next(); atoms != null; atoms = next()) {
            action.actionPerformed(atoms);
        }
    }

    public int size() {
        if(counter == null) {
            counter = new AtomsetCount();
        } else {
            counter.reset();
        }
        allAtoms(counter);
        return counter.callCount();
    }

    public int nBody() {
        return 3;
    }

    
    protected int stateUpDown; // 0 = up, 1 = middle, 2 = down, 3 = all done
    protected boolean doGoUp, doGoDown;
    protected AtomSet childList;
    protected IAtom targetAtom;
    protected int cursor;
    protected final AtomsetArray next;
    protected final IAtom[] atomArray;
    protected AtomsetCount counter;
    private static final long serialVersionUID = 1L;
}
