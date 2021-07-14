/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomSource;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveClusterAtomQ extends MCMoveAtom {

    protected final MCMoveAtom mcMove;
    
    public MCMoveClusterAtomQ(IRandom random, Space space, MCMoveAtom mcMove) {
        super(random, null, space);
        this.mcMove = mcMove;
    }

    public boolean doTrial() {
        boolean b = mcMove.doTrial();
        if (!b) return false;
        IAtomList leafAtoms = mcMove.getBox().getLeafList();
        int n = leafAtoms.size();
        for (int i=0; i<n; i++) {
            double u1 = random.nextDouble();
            double u2 = 2*Math.PI*random.nextDouble();
            double u3 = 2*Math.PI*random.nextDouble();
            double s1 = Math.sqrt(u1);
            double s2 = Math.sqrt(1-u1);
            Vector q = ((AtomOrientedQuaternion)leafAtoms.get(i)).getQuaternion();
            q.setX(0, s1*Math.sin(u2));
            q.setX(1, s1*Math.cos(u2));
            q.setX(2, s2*Math.sin(u3));
            q.setX(3, s2*Math.cos(u3));
        }
        return true;
    }

    public void setBox(Box p) {
        mcMove.setBox(p);
    }

    public double getChi(double temperature) {
        return mcMove.getChi(temperature);
    }

    public void rejectNotify() {
        mcMove.rejectNotify();
    }

    public void acceptNotify() {
        mcMove.acceptNotify();
    }

    public double energyChange() {
        return mcMove.energyChange();
    }

    public AtomIterator affectedAtoms() {
        return mcMove.affectedAtoms();
    }

    public AtomSource getAtomSource() {
        return mcMove.getAtomSource();
    }

    public void setAtomSource(AtomSource source) {
        mcMove.setAtomSource(source);
    }

    public Box getBox() {
        return mcMove.getBox();
    }

    public AtomIterator affectedAtoms(Box aBox) {
        return mcMove.affectedAtoms(aBox);
    }

    public double energyChange(Box aBox) {
        return mcMove.energyChange(aBox);
    }

    public String toString() {
        return mcMove.toString()+" with random quaternion rotation";
    }
    
    
}
