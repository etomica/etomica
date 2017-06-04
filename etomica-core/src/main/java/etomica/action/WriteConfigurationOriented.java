/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;

/**
 * WriteConfiguration subclass capable of handling oriented atoms.  Primary and
 * secondary directions are written on the same line as the atom's position.
 * 
 * @author Andrew Schultz
 */
public class WriteConfigurationOriented extends WriteConfiguration {

    public WriteConfigurationOriented(Space space) {
        super(space);
    }
    
    protected void writeAtom(FileWriter fileWriter, IAtom a) throws IOException {
        super.writeAtom(fileWriter, a);
        if (a instanceof IAtomOriented) {
            IOrientation o = ((IAtomOriented)a).getOrientation();
            Vector direction = o.getDirection();
            for (int i=0; i<direction.getD(); i++) {
                fileWriter.write(" "+direction.getX(i));
            }
            if (o instanceof IOrientationFull3D) {
                direction = ((IOrientationFull3D)o).getSecondaryDirection();
                for (int i=0; i<direction.getD(); i++) {
                    fileWriter.write(" "+direction.getX(i));
                }
            }
        }
    }
}
