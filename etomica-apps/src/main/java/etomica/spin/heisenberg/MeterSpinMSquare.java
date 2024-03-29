/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSource;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Undefined;


/**
 * returns the average of square of total dipole moment.
 *
 * @author Weisong Lin
 */
public class MeterSpinMSquare extends DataSourceScalar implements IDataSource {

    private final Vector sum;
    private final double dipoleMagnitude;
    private Box box;

    /**
     * @param space, box and dipoleMagnitude
     */

    public MeterSpinMSquare(Space space, Box box, double dipoleMagnitude) {
        super("Spin", Undefined.DIMENSION);
        sum = space.makeVector();
        this.box = box;
        this.dipoleMagnitude = dipoleMagnitude;
    }

    /**
     * @return <M^2> with M is the total dipole moment;
     */
    public double getDataAsScalar() {
        sum.E(0.0);
        for (IAtom a : box.getLeafList()) {
            sum.PE(((IAtomOriented)a).getOrientation().getDirection());
        }
        return sum.squared() * dipoleMagnitude * dipoleMagnitude;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }
}
