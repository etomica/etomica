/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Space;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

public class AtomLeafDynamic extends Atom implements IAtomKinetic {

    protected final Vector velocity;

    public AtomLeafDynamic(Space space, AtomType type) {
        super(space, type);
        velocity = space.makeVector();
    }

    public Vector getVelocity() {
        return velocity;
    }

    public void copyCoordinatesFrom(IAtom atom) {
        super.copyCoordinatesFrom(atom);
        velocity.E(((IAtomKinetic) atom).getVelocity());
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        super.saveState(fw);
        int D = velocity.getD();
        fw.write("" + velocity.getX(0));
        for (int i = 1; i < D; i++) fw.write(" " + velocity.getX(i));
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        super.restoreState(br);
        String[] bits = br.readLine().split(" ");
        int D = position.getD();
        for (int i = 0; i < D; i++) velocity.setX(i, Double.parseDouble(bits[i]));
    }
}
