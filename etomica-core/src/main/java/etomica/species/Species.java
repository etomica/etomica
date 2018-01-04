/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.config.IConformation;
import etomica.molecule.IMolecule;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Type for atom that is a group of other atoms, and for which its node is an instance of AtomTreeNodeGroup.
 *
 * @author andrew
 */
public abstract class Species implements ISpecies {

    protected final List<AtomType> atomTypes;
    protected int index;
    protected IConformation conformation;

    public Species() {
        this.atomTypes = new ArrayList<>();
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int newIndex) {
        index = newIndex;
    }

    public int getAtomTypeCount() {
        return this.atomTypes.size();
    }

    public AtomType getAtomType(int index) {
        return this.atomTypes.get(index);
    }

    @Override
    public List<AtomType> getAtomTypes() {
        return Collections.unmodifiableList(atomTypes);
    }

    public void initializeConformation(IMolecule molecule) {
        conformation.initializePositions(molecule.getChildList());
    }

    public void removeChildType(AtomType removedType) {
        boolean success = this.atomTypes.remove(removedType);
        if (!success) {
            throw new IllegalArgumentException("AtomType " + removedType + " is not my child!");
        }
    }

    public void addChildType(AtomType newChildType) {
        if (newChildType.getSpecies() != null) {
            throw new IllegalArgumentException(newChildType + " already has a parent");
        }
        newChildType.setSpecies(this);
        this.atomTypes.add(newChildType);
    }

    public IConformation getConformation() {
        return conformation;
    }

    public void setConformation(IConformation config) {
        conformation = config;
    }
}
