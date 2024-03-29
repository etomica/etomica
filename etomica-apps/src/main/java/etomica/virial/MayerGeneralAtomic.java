/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneralAtomic implements MayerFunction {

    /**
     * Constructor Mayer function using given potential.
     */
    public MayerGeneralAtomic(Space space, IPotential2 potential) {
        this.potential = potential;
        this.space = space;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        IAtom atom0 = pair.get(0).getChildList().get(0);
        IAtom atom1 = pair.get(1).getChildList().get(0);
        Vector dr = space.makeVector();
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        double x = -beta*potential.u(dr, atom0, atom1);
        if (Math.abs(x) < 0.01) {
            return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        return Math.exp(x) - 1;
    }

    public void setBox(Box newBox) {
    }

    protected final IPotential2 potential;
    protected final Space space;
}
