package etomica;

import etomica.atom.AtomSequencerFactory;
import etomica.atom.iterator.ApiInterspecies1A;
import etomica.atom.iterator.ApiInterspeciesAA;
import etomica.atom.iterator.ApiIntraspecies1A;
import etomica.atom.iterator.ApiIntraspeciesAA;
import etomica.atom.iterator.ApiMolecule;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;

public class IteratorFactorySimple extends IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
   
    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
    public AtomsetIteratorMolecule makeInterspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiInterspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiInterspeciesAA(species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiIntraspecies1A(species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAA(species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerFactory.SIMPLE;
    }
}