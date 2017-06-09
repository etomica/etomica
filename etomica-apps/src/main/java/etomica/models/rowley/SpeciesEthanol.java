/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;


import etomica.atom.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for ethanol with satellite site (Rowley et al 2006).
 */
public class SpeciesEthanol extends Species {

    public final static int indexO = 0;
    public final static int indexaC = 1;
    public final static int indexaH = 2; // ahType
    public final static int indexC = 3;
    public final static int indexH1a = 4; // hType
    public final static int indexH1b = 5; // hType
    public final static int indexH2a = 6; // hType
    public final static int indexH2b = 7; // hType
    public final static int indexH2c = 8; // hType
    public final static int indexX = 9;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType oType, acType, ahType, cType, hType, xType;
    public SpeciesEthanol(Space space, boolean pointCharges) {

        super();

        this.space = space;

        oType = new AtomType(Oxygen.INSTANCE);
        acType = new AtomType(Carbon.INSTANCE);
        ahType = new AtomType(Hydrogen.INSTANCE);
        cType = new AtomType(Carbon.INSTANCE);
        hType = new AtomType(Hydrogen.INSTANCE);
        xType = new AtomType(new ElementSimple("X", 1.0));

        addChildType(oType);
        addChildType(acType);
        addChildType(ahType);
        addChildType(cType);
        addChildType(hType);
        addChildType(xType);

        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationEthanol(space, pointCharges));
     }

     public IMolecule makeMolecule() {
         Molecule ethanol = new Molecule(this, 10);

         // The order in which the child atoms are added is important; it must match the order of site indices below.
         ethanol.addChildAtom(new Atom(space, oType));
         ethanol.addChildAtom(new Atom(space, acType));
         ethanol.addChildAtom(new Atom(space, ahType));
         ethanol.addChildAtom(new Atom(space, cType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, xType));
         conformation.initializePositions(ethanol.getChildList());
         return ethanol;
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

    public AtomType getCarbonType() {
         return cType;
     }

    public AtomType getHydrogenType() {
         return hType;
     }

    public AtomType getXType() {
         return xType;
     }

     public int getNumLeafAtoms() {
         return 10;
     }
}
