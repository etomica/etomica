/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;
import java.io.Serializable;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMoleculeList;
import etomica.atom.AtomsetArray;
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

    public void setBasis(IMoleculeList atoms) {
        if (atoms == null) {
            childList = null;
        }
        else {
            childList = atoms.getMolecule(0).getChildList();
        }
    }

    public int basisSize() {
        return 1;
    }

    public boolean haveTarget(IAtom target) {
        if (childList == null) {
            return false;
        }
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
        if (childList == null) {
            return;
        }

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
        if (childList == null) {
            return;
        }
        cursor = childList.getAtomCount();
    }

    public IAtomList next() {
        if (childList == null) {
            return null;
        }
        
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

    public int size() {
        if (childList == null) {
            return 0;
        }
        
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
    }

    public int nBody() {
        return 3;
    }

    
    protected int stateUpDown; // 0 = up, 1 = middle, 2 = down, 3 = all done
    protected boolean doGoUp, doGoDown;
    protected IAtomList childList;
    protected IAtom targetAtom;
    protected int cursor;
    protected final AtomsetArray next;
    protected final IAtom[] atomArray;
    private static final long serialVersionUID = 1L;
}
