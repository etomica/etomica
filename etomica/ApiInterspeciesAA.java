/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica;

/**
 * Gives pairs formed from the molecules of two different species in a phase.
 * Species are specified at construction and cannot be changed afterwards.
 */
public class ApiInterspeciesAA extends AtomsetIteratorAdapter implements
        AtomsetIteratorPhaseDependent {

    /**
     * @param species array of two different, non-null species
     */
    public ApiInterspeciesAA(Species[] species) {
        super(ApiBuilder.makeInterlistIterator());
        species0 = species[0];
        species1 = species[1];
        if(species0 == null || species1 == null) throw new NullPointerException("Constructor of ApiInterspeciesAA requires two non-null species");
        if(species0 == species1) throw new IllegalArgumentException("Constructor of ApiInterspeciesAA requires two different species");
        aiOuter = (AtomIteratorListSimple)((ApiInnerFixed)iterator).getOuterIterator();
        aiInner = (AtomIteratorListSimple)((ApiInnerFixed)iterator).getInnerIterator();
    }

    /** 
     * Configures iterator to return molecules from the set species in the given phase.
     */
    public void setPhase(Phase phase) {
        if(phase == null) {
            emptyList.clear();
            aiOuter.setList(emptyList);
        } else {
            aiOuter.setList(((AtomTreeNodeGroup)phase.getAgent(species0).node).childList);
            aiInner.setList(((AtomTreeNodeGroup)phase.getAgent(species1).node).childList);
        }
    }

    private final AtomIteratorListSimple aiInner, aiOuter;
    private final Species species0, species1;
    private final AtomList emptyList = new AtomList();
}
