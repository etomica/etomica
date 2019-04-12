/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;

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
public class SpeciesPolymerMono extends Species {

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomType leafAtomType;
    protected boolean isDynamic;
    private int f, l;

    /**
     * Constructs instance with a default element
     */
    public SpeciesPolymerMono(Simulation sim, Space _space, int f, int l) {
        this(_space, new ElementSimple(sim), f, l);
    }

    public SpeciesPolymerMono(Space _space, IElement element, int f, int l) {
        this(_space, new AtomType(element), f, l);
    }

    public SpeciesPolymerMono(Space space, AtomType leafAtomType, int f, int l) {
        super();
        this.space = space;
        this.leafAtomType = leafAtomType;
        addChildType(leafAtomType);
        this.f = f;
        this.l = l;
        setConformation(new ConformationStarPolymerGraft(space, f, l));
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
        Molecule group = new Molecule(this, f * l + 1);
        for (int i = 0; i < f * l + 1; i++) {
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
        return f * l + 1;
    }
}
