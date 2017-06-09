/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.atom.MoleculeOriented;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.SpeciesOriented;

public class SpeciesWater4POriented extends SpeciesOriented {

    public final static int indexH1 = 0;
    public final static int indexH2 = 1;
    public final static int indexO = 2;
    public final static int indexM = 3;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType oType, hType, mType;
    protected final boolean isDynamic;
    public SpeciesWater4POriented(Space space, boolean isDynamic) {
        super(space);
        this.isDynamic = isDynamic;
        this.space = space;
        hType = new AtomType(Hydrogen.INSTANCE);
        oType = new AtomType(Oxygen.INSTANCE);
        mType = new AtomType(new ElementSimple("M", 0.0));
        addChildType(hType);
        addChildType(oType);
        addChildType(mType);

        setConformation(new ConformationWaterTIP4P(space));
        init();
    }

    public IMolecule makeMolecule() {
        MoleculeOriented water = isDynamic ? new MoleculeOrientedDynamic(space, this, 4) :
                                             new MoleculeOriented(space, this, 4);
        water.addChildAtom(new Atom(space, hType));
        water.addChildAtom(new Atom(space, hType));
        water.addChildAtom(new Atom(space, oType));
        water.addChildAtom(new Atom(space, mType));
        conformation.initializePositions(water.getChildList());
        return water;
    }

    public AtomType getHydrogenType() {
        return hType;
    }

    public AtomType getOxygenType() {
        return oType;
    }

    public AtomType getMType() {
        return mType;
    }

    public int getNumLeafAtoms() {
        return 4;
    }
}
