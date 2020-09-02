/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.meta.annotations.IgnoreProperty;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Mass;

import java.util.StringJoiner;

/**
 * Identifies a set of atoms and defines properties of those atoms. Properties include indices used for tracking, mass
 * and element.
 */
public class AtomType {

    protected final IElement element;
    protected int index;
    protected ISpecies species;
    private final String name;

    public AtomType(IElement element) {
        this(element, element.getSymbol());
    }

    public AtomType(IElement element, String name) {
        this.element = element;
        this.name = name;
        this.index = -1;
    }

    /**
     * Convenience factory method for creating an AtomType with a one-off ElementSimple.
     *
     * @param symbol symbol for the element of this AtomType
     * @param mass   mass for the element of this AtomType
     * @return an AtomType with an ElementSimple
     */
    public static AtomType simple(String symbol, double mass) {
        return new AtomType(new ElementSimple(symbol, mass));
    }

    public static AtomType simple(String symbol) {
        return new AtomType(new ElementSimple(symbol));
    }

    public static AtomType simpleFromSim(Simulation sim) {
        return new AtomType(new ElementSimple(sim));
    }

    public static AtomType element(IElement element) {
        return new AtomType(element);
    }

    public static AtomType element(IElement element, String name) {
        return new AtomType(element, name);
    }

    /**
     * @return the index for this IAtomType, within the context of an Simulation.  The index is the IAtomType's position
     * in the list of atom types in the simulation.
     */
    public int getIndex() {
        return index;
    }

    /**
     *
     * @return the name of this AtomType. The name is by default the symbol of the element, but must be unique within
     * a Species.
     */
    public String getName() {
        return this.name;
    }

    /**
     * Informs the IAtomType what its index should be.  This should only be called by the species.
     *
     * @param newIndex the atom type's new index
     */
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    /**
     * @return the species that contains the atom type
     */
    @IgnoreProperty
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * Informs the atom type what species contains the atom types.  This should only be called by the species.
     *
     * @param newSpecies the atom type's new species
     */
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    /**
     * @return the value of the mass.
     */
    public final double getMass() {
        return element.getMass();
    }

    /**
     * @return the reciprocal of the mass, 1.0/mass
     */
    public final double rm() {
        return element.rm();
    }

    public final Dimension getMassDimension() {
        return Mass.DIMENSION;
    }

    /**
     * @return the element for this atom type
     */
    public final IElement getElement() {
        return element;
    }

    @Override
    public String toString() {
        return AtomType.class.getSimpleName() + "[" + this.name + "]";
    }
}
