/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;
import etomica.atom.IAtomType;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for 3-point water molecule.
 */
public class SpeciesWater3P extends Species {
    
    public SpeciesWater3P(Space space) {
        this(space, false);
    }
    
    public SpeciesWater3P(Space space, boolean isDynamic) {
        super();
        this.space = space;
        hType = new AtomTypeLeaf(Hydrogen.INSTANCE);
        oType = new AtomTypeLeaf(Oxygen.INSTANCE);
        addChildType(hType);
        addChildType(oType);
        this.isDynamic = isDynamic;

        setConformation(new ConformationWater3P(space));
    }
    
    public IMolecule makeMolecule() {
        Molecule water = new Molecule(this, 3);
        water.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));
        water.addChildAtom(isDynamic ? new AtomLeafDynamic(space, hType) : new Atom(space, hType));
        water.addChildAtom(isDynamic ? new AtomLeafDynamic(space, oType) : new Atom(space, oType));
        conformation.initializePositions(water.getChildList());
        return water;
    }
    
    public IAtomType getHydrogenType() {
        return hType;
    }
    
    public IAtomType getOxygenType() {
        return oType;
    }

    public int getNumLeafAtoms() {
        return 3;
    }
    
    public final static int indexH1 = 0;
    public final static int indexH2 = 1;
    public final static int indexO  = 2;

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomTypeLeaf oType, hType;
    protected final boolean isDynamic;
}
