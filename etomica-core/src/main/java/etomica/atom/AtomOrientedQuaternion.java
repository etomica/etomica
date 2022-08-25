/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.spaceNd.VectorND;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

public class AtomOrientedQuaternion extends Atom {

    protected final VectorND quaternion;

    public AtomOrientedQuaternion(Space space, AtomType type) {
        super(space, type);
        quaternion = new VectorND(4);
        quaternion.setX(0, 1);
    }

    public Vector getQuaternion() {
        return quaternion;
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        super.saveState(fw);
        fw.write("" + quaternion.getX(0));
        for (int i = 1; i < 4; i++) fw.write(" " + quaternion.getX(i));
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        super.restoreState(br);
        String[] bits = br.readLine().split(" ");
        for (int i = 0; i < 4; i++) position.setX(i, Double.parseDouble(bits[i]));
    }
}
