package etomica.atom.iterator;

import etomica.atom.ISpeciesAgent;
import etomica.box.Box;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of two different species in a box.
 * Species are specified at construction and cannot be changed afterwards.
 */

public class ApiInterspeciesAA extends AtomsetIteratorAdapter implements
        AtomsetIteratorBoxDependent {

    /**
     * Constructs iterator that provides iterates taken from the molecules of
     * two species. Given array is sorted in increasing order of species index. 
     * Then atoms in iterate pairs will be such that species of atom0 is
     * species[0] (having smaller species index), and species of atom1 is species[1]
     * (having larger species index).
     * 
     * @param species
     *            array of two different, non-null species
     * @throws NullPointerException
     *             if species or one of its elements is null
     * @throws IllegalArgumentException
     *             if species.length != 2 or if species[0] == species[1]
     */
    public ApiInterspeciesAA(Species[] species) {
        super(new ApiInterArrayList());
        apiInterList = (ApiInterArrayList) iterator;
        if(species.length != 2) {
            throw new IllegalArgumentException("Incorrect array length; must be 2 but length is "+species.length);
        }

        // we need to sort these.  we'll do that once we have the box
        species0 = species[0];
        species1 = species[1];
        if (species0 == null || species1 == null) {
            throw new NullPointerException(
                    "Constructor of ApiInterspeciesAA requires two non-null species");
        }
        if (species0 == species1) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspeciesAA requires two different species");
        }
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * box.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        ISpeciesAgent agent0 = box.getAgent(species0);
        ISpeciesAgent agent1 = box.getAgent(species1);
        if (agent0.getIndex() > agent1.getIndex()) {
            // species were out of order.  swap them
            Species tempSpecies = species0;
            species0 = species1;
            species1 = tempSpecies;
            ISpeciesAgent tempAgent = agent0;
            agent0 = agent1;
            agent1 = tempAgent;
        }
        apiInterList.setOuterList(agent0.getChildList());
        apiInterList.setInnerList(agent1.getChildList());
    }

    private static final long serialVersionUID = 1L;
    private final ApiInterArrayList apiInterList;
    private Species species0, species1;
}
