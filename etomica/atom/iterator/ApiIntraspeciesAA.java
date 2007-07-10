package etomica.atom.iterator;

import etomica.box.Box;
import etomica.species.Species;

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
    public ApiIntraspeciesAA(Species species) {
        this(new Species[] { species, species });
    }

    /**
     * @param species
     *            array of two non-null elements referencing the same species
     *            instance
     * @throws NullPointerException
     *             if species or one of its elements is null
     * @throws IllegalArgumentException
     *             if species is not a length-2 array or if its elements refer
     *             to different species instances
     */
    public ApiIntraspeciesAA(Species[] species) {
        super(new ApiIntraArrayList());
        if (species.length != 2)
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires two references to the same species instance");
        if (species == null || species[0] == null || species[1] == null)
            throw new NullPointerException(
                    "Constructor of ApiIntraspecies1A requires two non-null species references to the same instance");
        if (species[0] != species[1])
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires references to the same species instance");
        this.species = species[0];
        pairIterator = (ApiIntraArrayList) iterator;
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        pairIterator.setList(box.getMoleculeList(species));
    }

    private static final long serialVersionUID = 1L;
    private final ApiIntraArrayList pairIterator;
    private final Species species;
}
