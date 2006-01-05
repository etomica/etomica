package etomica.atom;

import etomica.atom.iterator.AtomIteratorArrayListSimple;


/**
 * Provides a way for AtomTreeNodeGroups to manage their ordinals.  They 
 * request new ordinals for newly inserted atoms and return ordinals for
 * removed atoms.  Returned ordinals are added to a reservoir, which is
 * used for new ordinals.  If an ordinal is returned and the reservoir is 
 * full, the ordinals from the reservoir are recycled are discarded.
 */
class OrdinalReservoir implements java.io.Serializable {

    public OrdinalReservoir(AtomTreeNodeGroup node) {
        atomNode = node;
        reservoirSize = 50;
        ordinalReservoir = new int[reservoirSize];
        maxOrdinal = 0;
        reservoirCount = 0;
    }

    /**
     * Returns a new ordinal.
     */
    public int requestNewOrdinal() {
        if (reservoirCount == 0) {
            return ++maxOrdinal;
        }
        reservoirCount--;
        return ordinalReservoir[reservoirCount];
    }
    
    /**
     * Indicates that this ordinal is no longer in use and can be discarded
     * or later reused for a new atom.
     */
    public void returnOrdinal(int ordinal) {
        if (ordinal == maxOrdinal) {
            return;
        }
        // add the index to the reservoir first
        ordinalReservoir[reservoirCount] = ordinal;
        reservoirCount++;
        // reservoir is full
        int oldReservoirCount = reservoirCount;
        AtomIteratorArrayListSimple childIterator = new AtomIteratorArrayListSimple(atomNode.childList);
        childIterator.reset();
        // loop over child atoms.  Any atoms whose ordinal exceeds the future
        // max is given a new one
        
        while (childIterator.hasNext()) {
            Atom a = childIterator.nextAtom();
            if (a.node.getOrdinal() > maxOrdinal-oldReservoirCount) {
                a.node.setOrdinal(requestNewOrdinal());
            }
        }
        maxOrdinal -= oldReservoirCount;
        if (reservoirCount != 0) {
            System.out.println("reservoir still has atoms:");
            for (int i=0; i<reservoirCount; i++) {
                System.out.print(ordinalReservoir[i]+" ");
            }
            throw new RuntimeException("I was fully expecting the reservoir to be empty!");
        }
    }
    
    public int getMaxGlobalOrdinal() {
        return maxOrdinal;
    }
    
    /**
     * Collapses the ordinals.  At the end, the maximum ordinal of 
     * a child atom will be maxOrdinal.  At the start, if ordinals 
     * in the reservoir are greater than the final maximum ordinal, they
     * are discarded.  Then Atoms with ordinals greater than the final 
     * maximum ordinal are given ordinals from the reservoir that are less
     * than the final maximum ordinal.  At the end maxOrdinal is decreased 
     * by the current previous number of ordinals in the reservoir.
     */
    private void collapseOrdinals() {
    }
    
    /**
     * Sets the size of the ordinal reservoir.
     */
    public void setReservoirSize(int size) {
        if (size < 0) {
            throw new IllegalArgumentException("Reservoir size must not be negative");
        }
        collapseOrdinals();
        // Set the actual reservoir size to one more because we collapse the
        // indices when it's full, not when it's full and we have another to add.
        reservoirSize = size+1;
        ordinalReservoir = new int[reservoirSize];
    }
    
    /**
     * Returns the size of the reservoir
     */
    public int getIndexReservoirSize() {
        return reservoirSize-1;
    }

    private int[] ordinalReservoir;
    private int reservoirSize;
    private int reservoirCount;
    private int maxOrdinal;
    private AtomTreeNodeGroup atomNode;
}
