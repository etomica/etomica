package etomica.nbr.cell;

import etomica.ApiIntergroup;
import etomica.AtomSequencerFactory;
import etomica.AtomsetIterator;
import etomica.IteratorFactory;
import etomica.Species;

public class IteratorFactoryCell extends IteratorFactory {
    
    public static final IteratorFactoryCell INSTANCE = new IteratorFactoryCell();
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomsetIterator makeInterSpeciesPairIterator(Species[] species) {
        return new ApiIntergroup();
    }
    public AtomsetIterator makeIntraSpeciesPairIterator(Species[] species) {
        return new ApiIntragroupCell(species);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}