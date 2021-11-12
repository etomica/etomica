/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.integrator.IntegratorBox;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Vector;

class DiameterHashVacancy implements DiameterHash {
    IntegratorBox integrator;
    double rc;
    double rc2;
    int nmax;
    NeighborIterator iter;

    public DiameterHashVacancy(IntegratorBox integrator, NeighborManager neighborManager, double rc) {
        this.integrator = integrator;
        this.rc = rc;
        rc2 = rc * rc;
        nmax = 12;
        iter = neighborManager.makeNeighborIterator();
    }

    public double getDiameter(IAtom a) {
        if (!integrator.getEventManager().firingEvent()) return 1.0;
        int[] n = {0};
        iter.iterAllNeighbors(a.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij) {
                if (rij.squared() < rc2) n[0]++;
            }
        });
        if (n[0] < nmax) return 1.0;
        if (n[0] == nmax) return 0.1;
        return 0.1;
    }
}
