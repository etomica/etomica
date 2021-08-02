/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSource;
import etomica.units.dimensions.Null;


/**
 * Meter that provides the x-component of the vector average of
 * spin values (which is represented by the atom's position
 * vector).
 *
 * @author David Kofke
 */
public class MeterSpinFasterer extends DataSourceScalar implements IDataSource {

    private final Box box;

    /**
     *
     */
    public MeterSpinFasterer(Box box) {
        super("Spin", Null.DIMENSION);
        this.box = box;
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.MeterScalar#getDataAsScalar(etomica.Box)
     */
    public double getDataAsScalar() {
        int count = 0;
        double sum = 0;
        for (IAtom atom : box.getLeafList()) {
            sum += ((IAtomOriented)atom).getOrientation().getDirection().getX(0);
            count++;
        }
        return sum / count;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
}
