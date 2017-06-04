/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;

/**
 * 
 * 
 * Species nitrogen molecule
 * 	with 4-points charges
 * 
 * @author Tai Boon Tan
 *
 */
public class SpeciesN2 extends Species {

    public final static int indexN1 = 0;
    public final static int indexN2 = 1;
    public final static int indexP1left = 2;
    public final static int indexP2left = 3;
    public final static int indexP1right = 4;
    public final static int indexP2right = 5;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType nType, pType;
    public SpeciesN2(Space space) {
        this(space, false);
    }
    public SpeciesN2(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        nType = new AtomType(Nitrogen.INSTANCE);
        pType = new AtomType(new ElementSimple("P", 1.0));
        addChildType(nType);
        addChildType(pType);

        setConformation(new ConformationNitrogen(space));
     }

     public IMolecule makeMolecule() {
         Molecule nitrogen = new Molecule(this, 6);
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogen.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));

         conformation.initializePositions(nitrogen.getChildList());
         return nitrogen;
     }

    public AtomType getNitrogenType() {
         return nType;
     }

    public AtomType getPType() {
         return pType;
     }

     public int getNumLeafAtoms() {
         return 6;
     }

     public void initializeConformation(IMolecule molecule, Vector v) {
         ((ConformationNitrogen)conformation).initializePositions(molecule.getChildList(), v);
     }
}
