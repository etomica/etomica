/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.traPPE;

import etomica.atom.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for methanol with satellite site (Rowley et al 2006).
 */
public class SpeciesMethanol extends Species {

    public final static int indexCH3 = 0;
    public final static int indexO = 1;
    public final static int indexH = 2;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType cH3Type, oType, hType;
    
    public SpeciesMethanol(Space space) {

        super();

        this.space = space;

        cH3Type = new AtomType(new ElementSimple("cH3", 1.0)); // diameter taken to be CH3-CH3 equilibrium LJ distance
        oType = new AtomType(Oxygen.INSTANCE); // diameter taken to be O-O equilibrium LJ distance
        hType = new AtomType(Hydrogen.INSTANCE); // H-H equilibrium distance is not applicable

        addChildType(cH3Type);
        addChildType(oType);
        addChildType(hType);

        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationMethanol(space));
     }

     public IMolecule makeMolecule() {
         Molecule methanol = new Molecule(this, 3);

         // The order in which the child atoms are added is important; it must match the site indices.
         methanol.addChildAtom(new Atom(space, cH3Type));
         methanol.addChildAtom(new Atom(space, oType));
         methanol.addChildAtom(new Atom(space, hType));

         conformation.initializePositions(methanol.getChildList());
         return methanol;
     }

    public AtomType getCH3Type() {
         return cH3Type;
     }

    public AtomType getOType() {
         return oType;
     }

    public AtomType getHType() {
         return hType;
     }

     public int getNumLeafAtoms() {
         return 3;
     }
}
