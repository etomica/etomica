/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveClusterAtomQ extends MCMoveBox {

    protected final IRandom random;
    protected final MCMoveBox mcMove;
    
    public MCMoveClusterAtomQ(IRandom random, Box box, MCMoveBox mcMove) {
        super();
        this.random = random;
        this.mcMove = mcMove;
        setBox(box);
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

    public Box getBox() {
        return mcMove.getBox();
    }

    public double energyChange(Box aBox) {
        return mcMove.energyChange(aBox);
    }

    public String toString() {
        return mcMove.toString()+" with random quaternion rotation";
    }
    
    
}
