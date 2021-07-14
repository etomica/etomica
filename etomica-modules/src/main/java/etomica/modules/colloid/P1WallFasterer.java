/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.colloid;

import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.potential.P1HardFieldGeneric;
import etomica.space.Vector;

public class P1WallFasterer extends P1HardFieldGeneric {

    protected final Box box;
    protected int chainLength;

    public P1WallFasterer(Box box, double[] collisionPositions, double[] energies, int chainLength) {
        super(1, collisionPositions, energies, true);
        this.box = box;
        this.chainLength = chainLength;
    }

    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }

    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        int idx = atom.getParentGroup().getIndex();
        int childIdx1 = idx % chainLength;
        if (childIdx1 > 0) {
            IAtomKinetic bondAtom = (IAtomKinetic) box.getMoleculeList(atom.getParentGroup().getType()).get(idx - 1).getChildList().get(0);
            double by = bondAtom.getPosition().getX(1);
            by += falseTime * bondAtom.getVelocity().getX(1);
            if (r.getX(fieldDimension) * by < 0) {
                // opposite sides of the box.  ignore this collision so that this atom can go back
                du[0] = 0;
                deltaP.E(0);
                return (r.getX(fieldDimension) > 0) ? 0 : collisionPositions.length - 2;
            }
        }
        return super.bump(atom, oldState, r, falseTime, deltaP, du);
    }
}
