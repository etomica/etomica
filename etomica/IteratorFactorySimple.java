package etomica;

public class IteratorFactorySimple extends IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomsetIterator makeInterSpeciesPairIterator(Species[] species) {
        return new ApiIntergroup();
    }
    public AtomsetIterator makeIntraSpeciesPairIterator(Species[] species) {
        return new ApiIntragroup();
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
}