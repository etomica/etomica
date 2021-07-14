/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.integrator;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.util.Debug;

public class ActionZeroMomentum implements IAction {

    protected final Box box;

    public ActionZeroMomentum(Box box) {
        this.box = box;
    }

    public void actionPerformed() {
        zeroMomenta(box);
    }

    public static void zeroMomenta(Box box) {
        Vector momentum = box.getSpace().makeVector();
        momentum.E(0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        if (nLeaf == 0) return;
        double totalMass = 0;
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtom a = leafList.get(iLeaf);
            double mass = a.getType().getMass();
            if (mass != Double.POSITIVE_INFINITY) {
                momentum.PEa1Tv1(mass, ((IAtomKinetic) a).getVelocity());
                totalMass += mass;
            }
        }
        if (totalMass == 0) return;
        momentum.TE(1.0 / totalMass);
        //momentum is now net velocity
        //set net momentum to 0
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            double rm = a.getType().rm();
            if (rm != 0 && rm != Double.POSITIVE_INFINITY) {
                a.getVelocity().ME(momentum);
            }
        }
        if (Debug.ON) {
            momentum.E(0);
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
                double mass = a.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass, a.getVelocity());
                }
            }
            momentum.TE(1.0 / totalMass);
            if (Math.sqrt(momentum.squared()) > 1.e-10) {
                System.out.println("Net momentum per leaf atom is " + momentum + " but I expected it to be 0");
            }
        }
        momentum.E(0);
    }
}
