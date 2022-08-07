/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

public class AtomOriented extends Atom implements
        IAtomOriented {

    protected final IOrientation iOrientation;

    public AtomOriented(Space space, AtomType type) {
        this(space, type, false);
    }

    public AtomOriented(Space space, AtomType type, boolean isAxisSymmetric) {
        super(space, type);
        if (space.D() == 3) {
            if (isAxisSymmetric) {
                iOrientation = new Orientation3D(space);
            }
            else {
                iOrientation = new OrientationFull3D(space);
            }
        }
        else {
            iOrientation = space.makeOrientation();
        }
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    public void copyCoordinatesFrom(IAtom atom) {
        super.copyCoordinatesFrom(atom);
        iOrientation.E(((IAtomOriented) atom).getOrientation());
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        super.saveState(fw);
        int D = position.getD();
        Vector direction = iOrientation.getDirection();
        fw.write("" + direction.getX(0));
        for (int i = 1; i < D; i++) fw.write(" " + direction.getX(i));
        if (iOrientation instanceof OrientationFull3D) {
            direction = ((OrientationFull3D) iOrientation).getSecondaryDirection();
            for (int i = 0; i < D; i++) fw.write(" " + direction.getX(i));
        }
        fw.write("\n");
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        super.restoreState(br);
        String[] bits = br.readLine().split(" ");
        int D = position.getD();
        Vector direction = iOrientation.getDirection();
        for (int i = 0; i < D; i++) direction.setX(i, Double.parseDouble(bits[1]));
        if (iOrientation instanceof OrientationFull3D) {
            Vector secondaryDirection = ((OrientationFull3D) iOrientation).getSecondaryDirection();
            for (int i = 0; i < D; i++) secondaryDirection.setX(i, Double.parseDouble(bits[3 + i]));
            ((OrientationFull3D) iOrientation).setDirections(direction, secondaryDirection);
        } else {
            iOrientation.setDirection(direction);
        }
    }
}
