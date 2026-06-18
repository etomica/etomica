/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.graphics.ColorScheme;
import etomica.integrator.IntegratorBox;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Vector;

import java.awt.*;

class ColorSchemeVacancy extends ColorScheme {
    IntegratorBox integrator;
    Vector dr;
    double rc;
    double rc2;
    int nmax;
    NeighborIterator iter;

    public ColorSchemeVacancy(IntegratorBox integrator, NeighborManager neighborManager, double rc) {
        this.integrator = integrator;
        dr = this.integrator.getBox().getSpace().makeVector();
        this.rc = rc;
        rc2 = rc * rc;
        nmax = 12;
        iter = neighborManager.makeNeighborIterator();
    }

    public Color getAtomColor(IAtom a) {
        if (!integrator.getEventManager().firingEvent()) return new Color(1.0f, 1.0f, 1.0f);

        final int[] n = {0};
        iter.iterAllNeighbors(a.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij, int nn) {
                if (rij.squared() < rc2) n[0]++;
            }
        });
        if (n[0] < nmax) return Color.RED;
        if (n[0] == nmax) return Color.BLUE;
        return Color.BLUE;
    }
}
