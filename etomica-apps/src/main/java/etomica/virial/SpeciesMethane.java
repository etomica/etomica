/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for methane with explicit hydrogen
 * Bond angle = 109.5
 * @shu
 * 01-27-2013
 */
public class SpeciesMethane extends Species {

    public final static int indexC = 0;
    public final static int indexH1 = 1;
    public final static int indexH2 = 2;
    public final static int indexH3 = 3;
    public final static int indexH4 = 4;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType cType, hType;

    public SpeciesMethane(Space space) {
        this(space, false);
    }
    public SpeciesMethane(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        cType = new AtomType(Carbon.INSTANCE);
        hType = new AtomType(new ElementSimple("cH3", 1.0));
        addChildType(cType);
        addChildType(hType);
        setConformation(new ConformationMethane(space));
     }

     public IMolecule makeMolecule() {
         Molecule methane = new Molecule(this, 5);

         // The order in which the child atoms are added is important; it must match the site indices.
         methane.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));

         methane.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));
         methane.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));
         methane.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));
         methane.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));

         conformation.initializePositions(methane.getChildList());

         return methane;

     }

    public AtomType getCType() {
         return cType;
     }

    public AtomType getHType() {
         return hType;
     }

    public int getNumLeafAtoms() {
         return 5;
     }
}
