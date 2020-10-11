/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.AtomTestCollective;
import etomica.atom.IAtom;
import etomica.data.IDataSource;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomTestLiquid implements AtomTestCollective {

    public AtomTestLiquid(Space space, IDataSource meterDeformation) {
        axis = space.makeVector();
        work = space.makeVector();
        meter = meterDeformation;
    }

    public void setCutoff(double newCutoff) {
        cutoff = newCutoff;
        cutoffSq = cutoff*cutoff;
    }
    
    public double getCutoff() {
        return cutoff;
    }

    public void resetTest() {
        double deformation = meter.getData().getValue(1);
        double factor = (1+deformation) / (1-deformation);
        axis.E(Math.pow(factor, -1.0/3.0));
        axis.setX(2,1.0/(axis.getX(0)*axis.getX(0)));
    }
    
    public boolean test(IAtom a) {
        Vector p = a.getPosition();
        work.E(p);
        work.DE(axis);
        double r2 = work.squared();
        return r2 < cutoffSq;
    }

    private static final long serialVersionUID = 1L;
    protected double cutoff, cutoffSq;
    protected final IDataSource meter;
    protected final Vector axis;
    protected final Vector work;
}
