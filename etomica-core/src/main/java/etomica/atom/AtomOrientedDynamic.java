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

public class AtomOrientedDynamic extends AtomLeafDynamic implements
        IAtomOrientedKinetic {

    protected final IOrientation iOrientation;
    protected final Vector angularVelocity;

    public AtomOrientedDynamic(Space space, AtomType type) {
        this(space, type, false);
    }

    public AtomOrientedDynamic(Space space, AtomType type, boolean isAxisSymmetric) {
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
        angularVelocity = space.makeVector();  //XXX wrong! angular velocity for 2D is a 1D vector
    }

    public Vector getAngularVelocity() {
        return angularVelocity;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    public void copyCoordinatesFrom(IAtom atom) {
        super.copyCoordinatesFrom(atom);
        iOrientation.E(((IAtomOriented) atom).getOrientation());
        angularVelocity.E(((IAtomOrientedKinetic) atom).getAngularVelocity());
    }

    public void saveState(Writer fw) throws IOException {
        super.saveState(fw);
        int D = position.getD();
        fw.write("" + angularVelocity.getX(0));
        for (int i = 1; i < D; i++) fw.write(" " + angularVelocity.getX(i));
        Vector direction = iOrientation.getDirection();
        for (int i = 0; i < D; i++) fw.write(" " + direction.getX(i));
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
        for (int i = 0; i < D; i++) angularVelocity.setX(i, Double.parseDouble(bits[i]));
        Vector direction = iOrientation.getDirection();
        for (int i = 0; i < D; i++) direction.setX(i, Double.parseDouble(bits[3 + i]));
        if (iOrientation instanceof OrientationFull3D) {
            Vector secondaryDirection = ((OrientationFull3D) iOrientation).getSecondaryDirection();
            for (int i = 0; i < D; i++) secondaryDirection.setX(i, Double.parseDouble(bits[6 + i]));
            ((OrientationFull3D) iOrientation).setDirections(direction, secondaryDirection);
        } else {
            iOrientation.setDirection(direction);
        }
    }
}
