/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IMolecule;
import etomica.space.Vector;
import etomica.atom.DipoleSource;
import etomica.atom.MoleculeOriented;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;

/**
 * Implementation of DipoleSource that can handle water molecules that are
 * MoleculeOriented.  The orientation is taken from the molecule rather than
 * calculated from scratch.
 * 
 * The dipole is assumed to point in the "secondary" direction held by the
 * orientation.  This is the correct direction for AtomWater3P.  The magnitude
 * of the dipole must be given to this class (via setDipoleStrength since that
 * depends on the potential in use.
 *
 * @author Andrew Schultz
 */
public class DipoleSourceWater implements DipoleSource {

    public DipoleSourceWater(Space space) {
        dipole = space.makeVector();
    }

    /**
     * Sets the strength (magnitude) of the dipole.
     */
    public void setDipoleStrength(double newDipoleStrength) {
        dipoleStrength = newDipoleStrength;
    }

    /**
     * Returns the dipole strength.
     */
    public double getDipoleStrength() {
        return dipoleStrength;
    }
    
    public Vector getDipole(IMolecule molecule) {
        // assume dipole points in the secondary orientation direction
        MoleculeOriented orientedMolecule = (MoleculeOriented)molecule;
        dipole.E(((IOrientationFull3D)orientedMolecule.getOrientation()).getSecondaryDirection());
        dipole.TE(dipoleStrength);
        return dipole;
    }

    protected double dipoleStrength;
    protected final Vector dipole;
}
