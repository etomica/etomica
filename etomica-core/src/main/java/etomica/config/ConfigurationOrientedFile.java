/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;

/**
 * ConfigurationFile subclass capable of handling oriented atoms.  When an
 * oriented atom is encountered, the primary and secondary directions of the
 * orientation are expected on the same line as the position.
 * 
 * @author Andrew Schultz
 */
public class ConfigurationOrientedFile extends ConfigurationFile {

    public ConfigurationOrientedFile(String aConfName, Space space) {
        super(aConfName);
        direction = space.makeVector();
        secondaryDirection = space.makeVector();
    }

    protected void setPosition(IAtom atom, String string) {
        super.setPosition(atom, string);
        if (atom instanceof IAtomOriented) {
            String[] coordStr = string.split(" +");
            for (int i=0; i<direction.getD(); i++) {
                direction.setX(i, Double.valueOf(coordStr[i+direction.getD()]).doubleValue());
            }
            IOrientation orientation = ((IAtomOriented)atom).getOrientation();
            if (orientation instanceof IOrientationFull3D) {
                for (int i=0; i<direction.getD(); i++) {
                    secondaryDirection.setX(i, Double.valueOf(coordStr[i+direction.getD()*2]).doubleValue());
                }
                ((IOrientationFull3D)orientation).setDirections(direction, secondaryDirection);
            }
            else {
                orientation.setDirection(direction);
            }
        }

    }
    
    private static final long serialVersionUID = 2L;
    protected final Vector direction, secondaryDirection;
}
