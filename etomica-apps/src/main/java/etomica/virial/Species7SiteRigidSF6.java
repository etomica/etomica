/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.chem.elements.Fluorine;
import etomica.chem.elements.Sulfur;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species SF6, 7 sites, LJ, rigid, no partial charge
 * Reference: Samios, Molecular force field investigation for sulfur hexafluoride: A computer simulation study
 * 
 * @author shu
 * 01-18-2013
 */
public class Species7SiteRigidSF6 {

    public final static int indexS = 0;
    public final static int indexF1 = 1;
    public final static int indexF2 = 2;
    public final static int indexF3 = 3;
    public final static int indexF4 = 4;
    public final static int indexF5 = 5;
    public final static int indexF6 = 6;

    public static SpeciesGeneral create(boolean isDynamic) {
        AtomType sType = new AtomType(Sulfur.INSTANCE);
        AtomType fType = new AtomType(Fluorine.INSTANCE);
        Space space = Space3D.getInstance();
        return new SpeciesBuilder(space)
                .setDynamic(isDynamic)
                .withConformation(new Conformation7SiteRigidSF6(space))
                .addCount(sType, 1)
                .addCount(fType, 6)
                .build();
    }
}
