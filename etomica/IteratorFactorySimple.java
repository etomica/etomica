package etomica;

public class IteratorFactorySimple extends IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomsetIteratorMolecule makeInterSpeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiInterspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiInterspeciesAA(species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraSpeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiIntraspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAA(species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
}