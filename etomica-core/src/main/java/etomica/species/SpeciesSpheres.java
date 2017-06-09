/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.api.IElement;
import etomica.atom.IMolecule;
import etomica.atom.*;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.config.IConformation;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */
public class SpeciesSpheres extends Species {

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType leafAtomType;
    protected boolean isDynamic;
    protected int atomsPerGroup;

    public SpeciesSpheres(Simulation sim, Space _space) {
        this(sim, _space, 1);
    }

    public SpeciesSpheres(Simulation sim, Space _space, int nA) {
        this(_space, nA, new ElementSimple(sim));
    }

    public SpeciesSpheres(Space _space, int nA, IElement leafElement) {
        this(nA, leafElement, new ConformationLinear(_space), _space);
    }
    
    public SpeciesSpheres(int nA, IElement leafElement,
    		              IConformation conformation, Space _space) {
        this(_space, nA, new AtomType(leafElement), conformation);
    }

    public SpeciesSpheres(Space _space, int nA, AtomType leafAtomType, IConformation conformation) {
        super();
        this.space = _space;
        addChildType(leafAtomType);
        setNumLeafAtoms(nA);
        setConformation(conformation);
        this.leafAtomType = leafAtomType;
    }

    public void setIsDynamic(boolean newIsDynamic) {
        isDynamic = newIsDynamic;
    }

    public boolean isDynamic() {
        return isDynamic;
    }

    public AtomType getLeafType() {
        return getAtomType(0);
    }

    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         Molecule group = new Molecule(this, atomsPerGroup);
         for(int i=0; i<atomsPerGroup; i++) {
             group.addChildAtom(makeLeafAtom());
         }
         conformation.initializePositions(group.getChildList());
         return group;
     }

    protected IAtom makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new Atom(space, leafAtomType);
     }

    public int getNumLeafAtoms() {
        return atomsPerGroup;
    }

    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     *
     * @param na The new number of atoms per group
     */
    public void setNumLeafAtoms(int na) {
        atomsPerGroup = na;
    }
}
