/*
 * History
 * Created on Nov 25, 2004 by kofke
 */
package etomica.nbr.cell;

/**
 * Interface indicating the object can reference a NeighborCell instance.
 * Expected use is in application to an atom sequencer (AtomSequencerCell) 
 * that holds a reference to the cell at its atom's current position.
 * @author kofke
 *
 */
public interface Cellular {

    /**
     * @return The current cell referenced by the instance.
     */
    public NeighborCell getCell();
}
