/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;

public class DiameterHashGlass extends DiameterHashByType {
    protected double fac = 1;
    protected boolean flipped = false;

    public void setFac(double fac) {
        this.fac = fac;
    }

    public double getFac() {
        return fac;
    }

    public double getDiameter(IAtom a) {
        return flipped ? 0 : getActualDiameter(a);
    }

    public void setFlipped(boolean flipped) {
        this.flipped = flipped;
    }

    public boolean getFlipped() {
        return flipped;
    }

    public double getActualDiameter(IAtom a) {
        return fac * super.getDiameter(a);
    }
}
