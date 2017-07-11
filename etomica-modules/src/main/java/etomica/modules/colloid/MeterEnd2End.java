/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Measures the end-to-end squared distance for the colloid system, averaged
 * across the grafted chains.
 * 
 * @author Andrew Schultz
 */
public class MeterEnd2End extends DataSourceScalar {

    protected Box box;
    protected int chainLength;
    protected Vector dr, drTot;
    
    public MeterEnd2End(Space space) {
        super("End2End", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}));
        dr = space.makeVector();
        drTot = space.makeVector();
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }

    public double getDataAsScalar() {
        IAtomList list = box.getLeafList();
        int nGraft = list.getAtomCount()/chainLength;
        int iAtom = 0;
        double sum = 0;
        for (int i=0; i<nGraft; i++) {
            drTot.E(0);
            dr.E(list.getAtom(0).getPosition());
            for (int j=1; j<chainLength; j++) {
                Vector p = list.getAtom(iAtom).getPosition();
                dr.ME(p);
                box.getBoundary().nearestImage(dr);
                drTot.PE(dr);
                
                dr.E(p);
                iAtom++;
            }
            sum += drTot.squared();
        }
        return sum/nGraft;
    }
}
