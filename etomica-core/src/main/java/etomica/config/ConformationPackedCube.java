/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Pack atoms sequentially on a cubic lattice.  Does not zero total momentum.
 *
 * @author David Kofke
 */

public class ConformationPackedCube implements IConformation {
    /**
     *
     * @param space species 2- or 3-dimensional space
     * @param bondLength distance between adjacent atoms on the lattice
     */
    public ConformationPackedCube(Space space, double bondLength) {
        this.space = space;
        this.bondLength = bondLength;
    }

    public void setBondLength(double b) {
        bondLength = b;
    }
    public double getBondLength() {return bondLength;}
    public Dimension getBondLengthDimension() {return Length.DIMENSION;}

    public void initializePositions(IAtomList atomList) {
        int n = atomList.size();
        if(n == 0) return;
        int count;
        if(space.D() == 2) {
            int gridSize = (int) Math.ceil(Math.sqrt(n));
            count = 0;
            for (int i = 0; i < gridSize && count < n; i += 2) {
                double x = bondLength * i;
                for (int j = 0; j < gridSize && count < n; j++) {
                    double y = bondLength * j;
                    atomList.get(count).getPosition().E(new double[]{x, y});
                    count++;
                }
                x = bondLength * (i + 1);
                for (int j = gridSize - 1; j >= 0 && count < n; j--) {
                    double y = bondLength * j;
                    atomList.get(count).getPosition().E(new double[]{x, y});
                    count++;
                }
            }
        } else {
            int gridSize = (int) Math.ceil(Math.cbrt(n));
            count = 0;
            for (int i = 0; i < gridSize && count < n; i++) {
                double x = bondLength * i;
                for (int j = 0; j < gridSize && count < n; j++) {
                    int nj = ( i%2 == 0 ) ? j : gridSize - 1 - j;
                    double y = bondLength * nj ;
                    for (int k = 0; k < gridSize && count < n; k++) {
                        double z = bondLength * (((i + nj) % 2 == 0) ? k : gridSize - 1 - k);
                         atomList.get(count).getPosition().E(new double[]{x, y, z});
                        count++;
                    }
                }
            }
        }
    }

    protected final Space space;
    protected double bondLength;
}
