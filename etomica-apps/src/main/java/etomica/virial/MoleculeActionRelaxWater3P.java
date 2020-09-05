/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * 
 */
package etomica.virial;

import etomica.action.MoleculeAction;
import etomica.atom.IAtom;
import etomica.models.water.SpeciesWater3P;
import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesGeneral;

public class MoleculeActionRelaxWater3P implements MoleculeAction {
    public MoleculeActionRelaxWater3P(Space space) {
        work = space.makeVector();
        cosAngle = Math.cos(109.5/180.0*Math.PI);
        sinAngle = Math.sin(109.5/180.0*Math.PI);
        distance = 1.0;
    }
    
    public void actionPerformed(IMolecule molecule) {
        ISpecies species = molecule.getType();
        IAtom O = molecule.getChildList().get(species.getAtomByTypeName("O"));
        IAtom H1 = molecule.getChildList().get(species.getAtomByTypeName("H", 1));
        IAtom H2 = molecule.getChildList().get(species.getAtomByTypeName("H", 2));
        // normalize OH1
        Vector p1 = H1.getPosition();
        p1.ME(O.getPosition());
        p1.TE(1/Math.sqrt(p1.squared()));
        Vector p2 = H2.getPosition();
        p2.ME(O.getPosition());
        p2.TE(1/Math.sqrt(p2.squared()));
        // move H2 to fix bond angle
        double d = p1.dot(p2);
        work.E(p2);
        work.PEa1Tv1(-d,p1);
        work.TE(1/Math.sqrt(work.squared()));
        p2.Ea1Tv1(sinAngle,work);
        p2.PEa1Tv1(cosAngle,p1);
        p2.TE(distance/Math.sqrt(p2.squared()));
        p1.TE(distance);
        p1.PE(O.getPosition());
        p2.PE(O.getPosition());
    }

    private final Vector work;
    private final double sinAngle, cosAngle, distance;
}
