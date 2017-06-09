/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.config.ConformationLinear;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Species in which molecules are each made of a single spherical atom.
 * Does not permit multiatomic molecules.  The advantage of this species
 * over the multiatomic version (used with 1 atom), is that one layer of
 * the atom hierarchy is eliminated in SpeciesSpheresMono.  Each atom is
 * the direct child of the species agent (i.e., each atom is at the "molecule"
 * level in the hierarchy, without an intervening AtomGroup).
 * 
 * @author David Kofke
 */
public class SpeciesSpheresMono extends Species {

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType leafAtomType;
    protected boolean isDynamic;

    /**
     * Constructs instance with a default element
     */
    public SpeciesSpheresMono(Simulation sim, Space _space) {
        this(_space, new ElementSimple(sim));
    }

    public SpeciesSpheresMono(Space _space, IElement element) {
        this(_space, new AtomType(element));
    }

    public SpeciesSpheresMono(Space space, AtomType leafAtomType) {
        super();
        this.space = space;
        this.leafAtomType = leafAtomType;
        addChildType(leafAtomType);
        setConformation(new ConformationLinear(space, 1));
    }

    public void setIsDynamic(boolean newIsDynamic) {
        isDynamic = newIsDynamic;
    }

    public boolean isDynamic() {
        return isDynamic;
    }

    public AtomType getLeafType() {
        return leafAtomType;
    }

    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         Molecule group = new Molecule(this, 1);
         group.addChildAtom(makeLeafAtom());
         return group;
     }

     protected IAtom makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new Atom(space, leafAtomType);
     }

     public int getNumLeafAtoms() {
         return 1;
     }
}
