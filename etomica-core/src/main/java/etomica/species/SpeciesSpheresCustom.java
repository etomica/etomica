/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.*;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.config.ConformationLinear;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.simulation.Simulation;
import etomica.space.Space;

import java.util.ArrayList;
import java.util.List;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * each of a specified type.
 * 
 * @author Andrew Schultz
 */

public class SpeciesSpheresCustom extends Species {

    protected Space space;
    protected boolean isDynamic;
    private int[] childAtomTypes;
    
    /**
     * Constructs instance with 0 components and total number of children
     * equal to 1.  The actual atom types must be set before use.
     */
    public SpeciesSpheresCustom(Simulation sim, Space _space) {
        this(sim,_space, 0);
    }
    
    /**
     * Constructs instance with the given number of atom types.  Generic atom
     * types are created.
     */
    public SpeciesSpheresCustom(Simulation sim, Space _space, int nTypes) {
        this(_space, makeElements(sim,nTypes));
    }
    
    /**
     * Constructs instance with the given elements.
     */
    public SpeciesSpheresCustom(Space _space, IElement[] leafElements) {
        this(_space, makeAtomTypes(leafElements));
    }

    /**
     * Constructs instance with the given atom types.
     */
    public SpeciesSpheresCustom(Space space, AtomType[] atomTypes) {
        super();
        this.space = space;
        for (int i=0; i<atomTypes.length; i++) {
            addChildType(atomTypes[i]);
        }
        setConformation(new ConformationLinear(space));
    }

    private static IElement[] makeElements(Simulation sim, int nTypes) {
        ElementSimple[] elements = new ElementSimple[nTypes];
        for (int i = 0; i < elements.length; i++) {
            elements[i] = new ElementSimple(sim);
        }
        return elements;
    }

    private static AtomType[] makeAtomTypes(IElement[] leafElements) {
        AtomType[] types = new AtomType[leafElements.length];
        for (int i = 0; i < types.length; i++) {
            types[i] = new AtomType(leafElements[i]);
        }
        return types;
    }

    public void setIsDynamic(boolean newIsDynamic) {
        isDynamic = newIsDynamic;
    }

    public boolean isDynamic() {
        return isDynamic;
    }

    /**
     * Constructs a new group containing a block of atoms for
     * each sub-type.
     */
    public IMolecule makeMolecule() {
        Molecule group = new Molecule(this, childAtomTypes.length);
        for (int childAtomTypeIndex : childAtomTypes) {
            group.addChildAtom(makeLeafAtom(atomTypes.get(childAtomTypeIndex)));
        }
        conformation.initializePositions(group.getChildList());
        return group;
    }

    private IAtom makeLeafAtom(AtomType leafType) {
        if (leafType instanceof AtomTypeOriented) {
            return isDynamic ? new AtomOrientedDynamic(space, leafType)
                             : new AtomOriented(space, leafType);
        }
        return isDynamic ? new AtomLeafDynamic(space, leafType)
                         : new Atom(space, leafType);
    }

    public void setChildAtomTypes(int[] childAtomTypes) {
        this.childAtomTypes = childAtomTypes;
    }

    public int getNumLeafAtoms() {
        return childAtomTypes.length;
    }
}
