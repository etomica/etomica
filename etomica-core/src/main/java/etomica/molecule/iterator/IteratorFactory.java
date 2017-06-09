/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.api.ISpecies;

/**
 * Class for construction of iterators of molecules.  Iterators are given for looping
 * over single molecules and pairs of molecules.  Different pair iterators are provided
 * for looping pairs from the same species, and pairs from different species.  
 * Straightforward iterators are given by IteratorFactorySimple, which generates pairs
 * from the childlists of species agents.  Iteration based on cell lists is performed by 
 * iterators given by IteratorFactoryCell, in the etomica.nbr.cell package.
 *
 * @author David Kofke
 */

//used by PotentialMaster.setSpecies

public class IteratorFactory implements java.io.Serializable {

    public static final IteratorFactory INSTANCE = new IteratorFactory();

    /**
     * Selects an appropriate iterator for the given species array.  If array contains
     * only one element, an atom iterator is returned. If array contains two elements,
     * an atom-pair iterator is returned, as given by the makeIntraSpeciesPairIterator
     * method if both elements of the array are equal, or as given by the
     * makeInterSpeciesPairIterator method if the array elements are different.
     * @param species array used to determine type of iterator to return
     * @return an appropriate iterator for looping over molecules of the given species
     */
    public MoleculesetIteratorPDT makeMoleculeIterator(ISpecies[] species) {
        if (species == null || species.length == 0 || species.length > 2
                || species[0] == null || species[species.length-1] == null) {
            throw new IllegalArgumentException("null or invalid number of species.  Must specify either 1 or 2 species instances.");
        }
        if (species.length==1) {
            return new MoleculeIteratorMolecule(species[0]);
        }
        if (species[0] == species[1]) {
            return makeIntraspeciesPairIterator(species[0]);
        }
        return makeInterspeciesPairIterator(species);
    }

    /**
     * creates a pair iterator which loops over all pairs in a neighbor list
     * between two groups
     * @return the pair iterator
     */
    public MoleculesetIteratorPDT makeInterspeciesPairIterator(ISpecies[] species) {
        MoleculesetIteratorPDT api1A = new MpiInterspecies1A(species);
        MoleculesetIteratorBoxDependent apiAA = new MpiInterspeciesAA(species);
        return new MpiMolecule(api1A, apiAA);
    }

    /**
     * creates a pair iterator which loops over all pairs in a neighbor list
     * within one group
     * @return the pair iterator
     */
    public MoleculesetIteratorPDT makeIntraspeciesPairIterator(ISpecies species) {
        MoleculesetIteratorPDT api1A = new MpiIntraspecies1A(species);
        MoleculesetIteratorBoxDependent apiAA = new MpiIntraspeciesAA(species);
        return new MpiMolecule(api1A, apiAA);
    }

    private static final long serialVersionUID = 1L;
}
