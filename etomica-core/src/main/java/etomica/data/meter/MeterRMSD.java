/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * Measured root mean squared displacement of all atoms in the box.  This meter
 * will properly track atoms across periodic boundaries so long as atoms do not
 * move more than L/2 between calls to the meter.
 * <p>
 * Created by andrew on 6/15/17.
 */
public class MeterRMSD extends DataSourceScalar {
    
    protected final Box box;
    protected final Vector[] originalPosition;
    protected final Vector[] lastPosition;
    protected final Vector dr, drTmp;
    
    public MeterRMSD(Box box, Space space) {
        super("RMSD", Length.DIMENSION);
        this.box = box;
        IAtomList atoms = box.getLeafList();
        originalPosition = new Vector[atoms.size()];
        lastPosition = new Vector[atoms.size()];
        for (int i = 0; i < atoms.size(); i++) {
            originalPosition[i] = space.makeVector();
            lastPosition[i] = space.makeVector();
            originalPosition[i].E(atoms.get(i).getPosition());
            lastPosition[i].E(originalPosition[i]);
        }
        dr = space.makeVector();
        drTmp = space.makeVector();
    }
    
    @Override
    public double getDataAsScalar() {
        double sum = 0;
        IAtomList atoms = box.getLeafList();
        Boundary boundary = box.getBoundary();
        for (int i = 0; i < atoms.size(); i++) {
            Vector p = atoms.get(i).getPosition();
            dr.Ev1Mv2(p, lastPosition[i]);
            drTmp.E(dr);
            boundary.nearestImage(drTmp);
            if (!dr.equals(drTmp)) {
                originalPosition[i].PE(dr);
                originalPosition[i].ME(drTmp);
            }
            sum += p.Mv1Squared(originalPosition[i]);
            lastPosition[i].E(p);
        }
        return Math.sqrt(sum / atoms.size());
    }
}
