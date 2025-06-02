/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Assumes cubic box and 3D, must wrap a spherical potential
 */
public class P2LatticeSum implements IPotential2 {

    protected final IPotential2 p2;
    protected Boundary boundary;

    // wrapped potential p2 must be spherical
    public P2LatticeSum(IPotential2 p2) {
        this.p2 = p2;
    }

    public void setBoundary(Boundary boundary) {
        this.boundary = boundary;
        Vector L = boundary.getBoxSize();
        if (L.getD() != 3) {
            throw new RuntimeException("Must be a 3D system");
        }
        if (L.getX(0) != L.getX(1) || L.getX(1) != L.getX(2)) {
            throw new RuntimeException("Box must be cubic");
        }
    }

    public double getRange() {
        return p2.getRange();
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        Vector L = boundary.getBoxSize();
        int nImages = (int)Math.floor(p2.getRange() / L.getX(0) + 0.5);
        Vector idr12 = new Vector3D();
        double sum = 0;
        for (int ix=-nImages; ix<=nImages; ix++) {
            idr12.setX(0, dr12.getX(0) + ix * L.getX(0));
            for (int iy = -nImages; iy <= nImages; iy++) {
                idr12.setX(1, dr12.getX(1) + iy * L.getX(1));
                for (int iz = -nImages; iz <= nImages; iz++) {
                    idr12.setX(2, dr12.getX(2) + iz * L.getX(2));
                    sum += p2.u(idr12.squared());
                }
            }
        }
        return sum;
    }
}
