/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * measures the maximum amount by which by which the HS diameter could increase
 * without causing overlap
 * 
 * @author Andrew Schultz
 */
public class MeterMaxExpansion extends DataSourceScalar {

    protected final Vector dr;
    protected final NeighborIterator neighborIterator;
    protected final Box box;

    public MeterMaxExpansion(Box box, NeighborIterator iterator) {
        super("displacement", Length.DIMENSION);
        this.neighborIterator = iterator;
        this.box = box;
        dr = box.getSpace().makeVector();
    }
    
    public double getDataAsScalar() {
        Boundary boundary = box.getBoundary();
        IAtomList leafList = box.getLeafList();
        double[] min = {1e10};
        for (int i = 0; i<leafList.size(); i++) {
            IAtom atomi = leafList.get(i);
            neighborIterator.iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    dr.Ev1Mv2(atomi.getPosition(), jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    min[0] = Math.min(min[0], r2);
                }
            });
        }
        return Math.sqrt(min[0])-1;
    }
}
