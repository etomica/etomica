package etomica;

import etomica.AtomSequencer.Factory;

public class IteratorFactorySimple extends IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public Factory interactionAtomSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
    public Factory interactionMoleculeSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
    public AtomsetIterator makeInterSpeciesPairIterator(Species[] species) {
        return new ApiIntergroup();
    }
    public AtomsetIterator makeIntraSpeciesPairIterator(Species[] species) {
        return new ApiIntragroup();
    }
    public Factory moleculeSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
}