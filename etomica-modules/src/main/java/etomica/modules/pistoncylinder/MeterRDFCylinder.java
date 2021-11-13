/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.box.Box;
import etomica.data.IData;
import etomica.data.meter.MeterRDF;
import etomica.space.Vector;

/**
 * MeterRDF sublcass that properly calculates the RDF for the Piston/Cylinder
 * apparatus.  Including all pairs would yield a very non-bulk RDF since atoms
 * near the wall have open space next tho them.  To account for this, if
 * "padding" is the collision radius + RDF cutoff, then the RDF calculation
 * includes all pairs where both are a distance "padding" away from all walls
 * and half the pairs where only one is a distance "padding" away from all
 * walls.  The RDF is normalized to account for this.
 *
 * @author Andrew Schultz
 */
public class MeterRDFCylinder extends MeterRDF {

    public MeterRDFCylinder(Box box, P1HardMovingBoundary pistonPotential) {
        super(box.getSpace());
        setBox(box);
        this.pistonPotential = pistonPotential;
        reset();
    }

    public void reset() {
        super.reset();
        // make a new iterator with a new filter.  xMax might have changed
        ApiFilteredCylinder.AtomFilterInCylinder filter = new ApiFilteredCylinder.AtomFilterInCylinder(box.getBoundary(), pistonPotential, xDataSource.getXMax());
        iterator = new ApiFilteredCylinder(filter);
    }

    public IData getData() {
        super.getData();

        // renormalize the RDF to account for the excluded pairs
        double pistonRatio = 1;
        Vector dimensions = box.getBoundary().getBoxSize();
        double radius = pistonPotential.getCollisionRadius();
        for (int i = 0; i < space.D(); i++) {
            if (i == 1) {
                double ySize = dimensions.getX(1) * 0.5;
                if (dimensions.getD() == 2) {
                    ySize -= pistonPotential.getWallPosition();
                } else {
                    ySize += pistonPotential.getWallPosition();
                }
                pistonRatio *= (ySize - 2 * (xMax + radius)) / (ySize - 2 * radius);
                if (pistonRatio < 0) {
                    pistonRatio = Double.NaN;
                }
            } else {
                pistonRatio *= (dimensions.getX(i) - 2 * (xMax + radius)) / (dimensions.getX(i) - 2 * radius);
            }
        }
        data.TE(1 / pistonRatio);
        return data;
    }

    protected etomica.modules.pistoncylinder.P1HardMovingBoundary pistonPotential;
}
