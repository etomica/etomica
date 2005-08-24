package etomica.atom.iterator;

import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;
import etomica.species.Species;

public class IteratorFactorySimple extends IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomLinker.FACTORY;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomLinker.FACTORY;
    }
    public AtomsetIteratorMolecule makeInterspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiInterspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiInterspeciesAA(species);
        ApiSpecies11 api11 = new ApiSpecies11(species);
        return new ApiMolecule(api11, api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiIntraspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAA(species);
        ApiSpecies11 api11 = new ApiSpecies11(species);
        return new ApiMolecule(api11, api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomLinker.FACTORY;
    }
}
