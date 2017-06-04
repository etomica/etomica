/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.IEtomicaDataSource;
import etomica.potential.P1Harmonic;
import etomica.units.Energy;


/**
 * Computes the free-enery difference through free-energy perturbation.
 * Written specifically for harmonic 1-body potential, but wouldn't be hard
 * to modify for more general cases.
 *
 * @author David Kofke
 *
 */
public class MeterFreeEnergy extends DataSourceScalar implements IEtomicaDataSource {
    
    public MeterFreeEnergy(P1Harmonic reference, P1Harmonic target) {
        super("Free energy", Energy.DIMENSION);
        this.reference = reference;
        this.target = target;
    }
    
    public double getDataAsScalar() {
        iterator.reset();
        double sum = 0.0;
        for (IAtomList a = iterator.next(); a != null; a = iterator.next()) {
            sum += target.energy(a) - reference.energy(a);
        }
        return Math.exp(-sum);
    }
    
    public void setBox(Box box) {
        iterator.setBox(box);
        this.box = box;
    }
    
    public Box getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    Box box;
    P1Harmonic reference, target;
}
