/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSource;
import etomica.potential.P1Harmonic;
import etomica.units.dimensions.Energy;


/**
 * Computes the free-enery difference through free-energy perturbation.
 * Written specifically for harmonic 1-body potential, but wouldn't be hard
 * to modify for more general cases.
 *
 * @author David Kofke
 *
 */
public class MeterFreeEnergy extends DataSourceScalar implements IDataSource {
    
    public MeterFreeEnergy(P1Harmonic reference, P1Harmonic target) {
        super("Free energy", Energy.DIMENSION);
        this.reference = reference;
        this.target = target;
    }
    
    public double getDataAsScalar() {
        double sum = 0.0;
        for (IAtom a : box.getLeafList()) {
            sum += target.u(a) - reference.u(a);
        }
        return Math.exp(-sum);
    }
    
    public void setBox(Box box) {
        this.box = box;
    }
    
    public Box getBox() {
        return box;
    }

    Box box;
    P1Harmonic reference, target;
}
