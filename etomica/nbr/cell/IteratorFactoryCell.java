package etomica.nbr.cell;

import etomica.AtomSequencerFactory;
import etomica.AtomsetIterator;
import etomica.AtomsetIteratorPhaseDependent;
import etomica.IteratorFactory;
import etomica.Species;

public class IteratorFactoryCell extends IteratorFactory {

    private final int D;
    
    /**
     * @param D dimension of the space 
     */
    public IteratorFactoryCell(int D) {
        this.D = D;
    }
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomsetIterator makeInterSpeciesPairIterator(Species[] species) {
        AtomsetIteratorPhaseDependent api1A = null;//new ApiIntergroupCell1A(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntergroupCellAA(D, species);
        return new ApiCellAdapter(api1A, apiAA);
    }
    public AtomsetIterator makeIntraSpeciesPairIterator(Species[] species) {
        AtomsetIteratorPhaseDependent api1A = null;//new ApiIntragroupCell1A(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntragroupCellAA(D, species);
        return new ApiCellAdapter(api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}