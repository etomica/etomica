/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.box.storage.DoubleStorage;
import etomica.space.Space;

public class AtomHydrogen extends AtomOriented {    
    protected DoubleStorage.DoubleWrapper bondLength;

    public AtomHydrogen(Space space, AtomTypeOriented atype, DoubleStorage.DoubleWrapper bondLength) {
        super(space, atype);        
//        bondLength = BohrRadius.UNIT.toSim(1.401065676);
        this.bondLength = bondLength;//BohrRadius.UNIT.toSim(1.448736);
    }

    public static double getAvgBondLength(double x) {
        double a0 = 0.766099;//0.0766111;
        double a1 = 6.49168e-06;//6.84704E-07;
        double a2 = 1.55595e-09;//-3.89889E-11;
        double a3 = -1.27347e-12;//6.73971E-14;
        double a4 = 5.2509e-16;
        double y = a0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x;
//        System.out.println(x+" "+y);
//        System.exit(1);
        return y;
    }

    public double getBondLength() {
        return bondLength.get();
    }

    public void setBondLength(double x) {
        bondLength.set(x);
    }

    @Override
    public void copyFrom(IAtom atom) {
        super.copyFrom(atom);
        this.setBondLength(((AtomHydrogen) atom).getBondLength());
    }
}
