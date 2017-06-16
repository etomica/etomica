/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationChainZigZag2;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresHetero;

public class SpeciesAlkane extends SpeciesSpheresHetero {

    // parameters for TraPPE
    protected static final double nominalBondL = 1.54;
    protected static final double nominalBondTheta = Math.PI * 114 / 180;
    private static final long serialVersionUID = 1L;

    public SpeciesAlkane(Space _space, int numCarbons) {
        this(_space,numCarbons, new ElementSimple("CH3", 15), new ElementSimple("CH2", 14));
    }

    public SpeciesAlkane(Space _space, int numCarbons, ElementSimple CH3element, ElementSimple CH2element) {
    	super(_space, makeAtomTypes(new Element[]{CH3element, CH2element}));
        setTotalChildren(numCarbons);
        setConformationParameters(nominalBondL, nominalBondTheta);
    }
    
    public void setConformationParameters(double bondL, double bondTheta) {
        Vector vector1 = space.makeVector();
        vector1.setX(0, bondL);
        Vector vector2 = space.makeVector();
        vector2.setX(0, -bondL*Math.cos(bondTheta));
        vector2.setX(1, bondL*Math.sin(bondTheta));
        conformation = new ConformationChainZigZag2(space, vector1, vector2);
    }

    public IMolecule makeMolecule() {
        Molecule group = new Molecule(this, totalChildCount);
        //make straight alkane CH3-CH2-...-CH2-CH3
        group.addChildAtom(makeLeafAtom(childTypes[0]));
        for(int j = 0; j < childCount[1]; j++) {
            group.addChildAtom(makeLeafAtom(childTypes[1]));//CH2
        }
        if (childCount[0] > 1) {
            group.addChildAtom(makeLeafAtom(childTypes[0]));//CH3
        }
        conformation.initializePositions(group.getChildList());
        return group;
    }

    public AtomType getCH3Type() {
        return childTypes[0];
    }

    public AtomType getCH2Type() {
        return childTypes[1];
    }

    public void setTotalChildren(int newTotalChildren) {
        if (newTotalChildren > 1) {
            childCount[0] = 2;// CH3
            childCount[1] = newTotalChildren - 2;//CH2
        }
        else {
            childCount[0] = 1;// CH4???????????????????????
            childCount[1] = 0;
        }
    }
}
