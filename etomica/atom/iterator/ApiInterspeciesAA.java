package etomica.atom.iterator;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of two different species in a phase.
 * Species are specified at construction and cannot be changed afterwards.
 */

public class ApiInterspeciesAA extends AtomPairIteratorAdapter implements
        AtomsetIteratorPhaseDependent {

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

        // we need to sort these.  we'll do that once we have the phase
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
     * phase.  No iterates are given if phase is null.
     */
    public void setPhase(Phase phase) {
        if (phase == null) {
            emptyList.clear();
            apiInterList.setOuterList(emptyList);
        } else {
            AtomTreeNodeGroup agentNode0 = (AtomTreeNodeGroup) phase.getAgent(species0).getNode();
            AtomTreeNodeGroup agentNode1 = (AtomTreeNodeGroup) phase.getAgent(species1).getNode();
            if (agentNode0.getIndex() > agentNode1.getIndex()) {
                // species were out of order.  swap them
                Species tempSpecies = species0;
                species0 = species1;
                species1 = tempSpecies;
                AtomTreeNodeGroup tempNode = agentNode0;
                agentNode0 = agentNode1;
                agentNode1 = tempNode;
            }
            apiInterList.setOuterList(agentNode0.getChildList());
            apiInterList.setInnerList(agentNode1.getChildList());
        }
    }

    private static final long serialVersionUID = 1L;
    private final ApiInterArrayList apiInterList;
    private Species species0, species1;
    private final AtomArrayList emptyList = new AtomArrayList();
}
