package etomica;

import etomica.AtomSequencer.Factory;

public class IteratorFactorySimple implements IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public Factory interactionAtomSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
    public Factory interactionMoleculeSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
    public AtomsetIterator makeInterSpeciesPairIterator() {
        return new ApiIntergroup();
    }
    public AtomsetIterator makeIntraSpeciesPairIterator() {
        return new ApiIntragroup();
    }
    public Factory moleculeSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
}