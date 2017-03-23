/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;
import etomica.api.IAtom;
import etomica.api.IAtomType;
import etomica.api.IElement;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomTypeOriented;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.simulation.Simulation;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * each of a specified type.
 * 
 * @author Andrew Schultz
 */

public class SpeciesSpheresCustom extends Species {

    /**
     * Constructs instance with 0 components and total number of children 
     * equal to 1.  The actual atom types must be set before use.
     */
    public SpeciesSpheresCustom(Simulation sim, ISpace _space) {
        this(sim,_space, 0);
    }
    
    /**
     * Constructs instance with the given number of atom types.  Generic atom
     * types are created.
     */
    public SpeciesSpheresCustom(Simulation sim, ISpace _space, int nTypes) {
        this(_space, makeElements(sim,nTypes));
    }
    
    private static IElement[] makeElements(Simulation sim, int nTypes) {
        ElementSimple[] elements = new ElementSimple[nTypes];
        for (int i=0; i<elements.length; i++) {
            elements[i] = new ElementSimple(sim);
        }
        return elements;
    }
    
    /**
     * Constructs instance with the given elements.
     */
    public SpeciesSpheresCustom(ISpace _space, IElement[] leafElements) {
        this(_space, makeAtomTypes(leafElements));
    }
    
    protected static final AtomTypeLeaf[] makeAtomTypes(IElement[] leafElements) {
        AtomTypeLeaf[] types = new AtomTypeLeaf[leafElements.length];
        for (int i=0; i<types.length; i++) {
            types[i] = new AtomTypeLeaf(leafElements[i]);
        }
        return types;
    }
    
    /**
     * Constructs instance with the given atom types.
     */
    public SpeciesSpheresCustom(ISpace space, IAtomType[] atomTypes) {
        super();
        this.space = space;
        for (int i=0; i<atomTypes.length; i++) {
            addChildType(atomTypes[i]);
        }
        setConformation(new ConformationLinear(space));
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
        Molecule group = new Molecule(this, atomTypes.length);
        for (int i=0; i<atomTypes.length; i++) {
            group.addChildAtom(makeLeafAtom(childTypes[atomTypes[i]]));
        }
        conformation.initializePositions(group.getChildList());
        return group;
    }

    protected IAtom makeLeafAtom(IAtomType leafType) {
        if (leafType instanceof IAtomTypeOriented) {
            return isDynamic ? new AtomOrientedDynamic(space, leafType)
                             : new AtomOriented(space, leafType);
        }
        return isDynamic ? new AtomLeafDynamic(space, leafType)
                         : new Atom(space, leafType);
    }

    public int getNumComponents() {
        return childTypes.length;
    }
    
    /**
     * Sets the factories that make the child atoms of this factory's atom.
     * If the number of factories changes, the number fractions are set so 
     * that there is an equal amount of each child.  The caller is responsible
     * for ensuring that the AtomTypes for the child factories are children
     * of this AtomFactory's AtomType.
     *
     * @throws IllegalArgumentException
     *             if newChildFactory is an empty array
     */
    public void setChildTypes(IAtomType[] newchildTypes) {
        for (int i=0; i<childTypes.length; i++) {
            removeChildType(childTypes[i]);
        }
        for (int i=0; i<childTypes.length; i++) {
            addChildType(newchildTypes[i]);
        }
    }
    
    public void setAtomTypes(int[] atomTypes) {
        this.atomTypes = atomTypes;
    }

    public int getNumLeafAtoms() {
        return atomTypes.length;
    }

    private static final long serialVersionUID = 1L;
    protected ISpace space;
    protected boolean isDynamic;
    protected int[] atomTypes;
}
