/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;


import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for methanol with satellite site (Rowley et al 2006).
 */
public class SpeciesMethanol extends Species {

    public final static int indexO = 0;
    public final static int indexaC = 1;
    public final static int indexaH = 2; // ahType
    public final static int indexH1 = 3; // hType
    public final static int indexH2a = 4; // hType
    public final static int indexH2b = 5; // hType
    public final static int indexX = 6;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType oType, acType, ahType, hType, xType;
    public SpeciesMethanol(Space space, boolean pointCharges) {

        super();

        this.space = space;

        oType = new AtomType(Oxygen.INSTANCE);
        acType = new AtomType(Carbon.INSTANCE);
        ahType = new AtomType(Hydrogen.INSTANCE);
        hType = new AtomType(Hydrogen.INSTANCE);
        xType = new AtomType(new ElementSimple("X", 1.0));

        addChildType(oType);
        addChildType(acType);
        addChildType(ahType);
        addChildType(hType);
        addChildType(xType);

        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationMethanol(space, pointCharges));
     }

     public IMolecule makeMolecule() {
         Molecule methanol = new Molecule(this, 7);

         // The order in which the child atoms are added is important; it must match the site indices.
         methanol.addChildAtom(new Atom(space, oType));
         methanol.addChildAtom(new Atom(space, acType));
         methanol.addChildAtom(new Atom(space, ahType));
         methanol.addChildAtom(new Atom(space, hType));
         methanol.addChildAtom(new Atom(space, hType));
         methanol.addChildAtom(new Atom(space, hType));
         methanol.addChildAtom(new Atom(space, xType));
         conformation.initializePositions(methanol.getChildList());
         return methanol;
     }

    public AtomType getOxygenType() {
         return oType;
     }

    public AtomType getAlphaCarbonType() {
         return acType;
     }

    public AtomType getAlphaHydrogenType() {
         return ahType;
     }

    public AtomType getHydrogenType() {
         return hType;
     }

    public AtomType getXType() {
         return xType;
     }

     public int getNumLeafAtoms() {
         return 7;
     }
}
