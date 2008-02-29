package etomica.atom.iterator;

import etomica.api.IBox;
import etomica.species.ISpecies;

/**
 * Gives pairs formed from the molecules of a single species in a box. Species
 * is specified at construction and cannot be changed afterwards.
 */
public class ApiIntraspeciesAA extends AtomsetIteratorAdapter implements
        AtomsetIteratorBoxDependent {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public ApiIntraspeciesAA(ISpecies species) {
        super(new ApiIntraArrayList());
        if (species == null) {
            throw new NullPointerException("Constructor of ApiIntraspecies1A requires a non-null species");
        }
        this.species = species;
        pairIterator = (ApiIntraArrayList) iterator;
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(IBox box) {
        pairIterator.setList(box.getMoleculeList(species));
    }

    private static final long serialVersionUID = 1L;
    private final ApiIntraArrayList pairIterator;
    private final ISpecies species;
}
