/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.AlkaneEH;


import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for TraPPE-Explicit hydrogen,Siepmann
 *
 * @author shu 01-30-2013
 */
public class SpeciesAlkaneEH extends Species {

    public final int numCarbons;
    public final int numH;
    public final int numCH2;
    protected final Space space;
    private final int totalAtoms;
    private final AtomType c_3Type, c_2Type, hType;
    public int indexC_3_1;//CH3 on the left
    public int indexC_3_2;//CH3 on the right
    protected boolean isDynamic;

    public SpeciesAlkaneEH(Space space, int numCarbons) {
        super();
        this.space = space;
        this.numCarbons = numCarbons;
        indexC_3_1 = 0;
        indexC_3_2 = numCarbons - 1;

        numCH2 = numCarbons - 2;
        numH = numCarbons * 2 + 2;
        totalAtoms = numCarbons * 3 + 2;

        c_3Type = AtomType.simple("c_3", 1.0);
        c_2Type = AtomType.simple("c_2", 1.0);
        hType = AtomType.simple("h", 1.0);
        addChildType(c_3Type);
        addChildType(c_2Type);
        addChildType(hType);

        setConformation(new ConformationAlkaneEH(space, this));
    }

    public IMolecule makeMolecule() {

        Molecule molecule = new Molecule(this, totalAtoms);

        molecule.addChildAtom(makeLeafAtom(c_3Type)); // the left C of CH3, global index [0]

        for (int j = 0; j < numCH2; j++) {//CH2, global index: [1]~[n-2]
            molecule.addChildAtom(makeLeafAtom(c_2Type));
        }

        molecule.addChildAtom(makeLeafAtom(c_3Type)); // the right C of CH3,global index [n-1]

        for (int i = 0; i < numH; i++) {//H
            molecule.addChildAtom(makeLeafAtom(hType));
        }
        conformation.initializePositions(molecule.getChildList());
        return molecule;
    }

    public AtomType getC_3Type() {
        return this.c_3Type;
    }

    public AtomType getC_2Type() {
        return this.c_2Type;
    }

    public AtomType getHType() {
        return this.hType;
    }

    public int getNumLeafAtoms() {
        return totalAtoms;
    }

    protected IAtom makeLeafAtom(AtomType leafType) {
        return isDynamic ? new AtomLeafDynamic(space, leafType) : new Atom(space, leafType);
    }

}
