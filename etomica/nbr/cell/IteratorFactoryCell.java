package etomica.nbr.cell;

import etomica.ApiInnerVariable;
import etomica.ApiIntergroup;
import etomica.ApiIntragroup;
import etomica.AtomIteratorBasis;
import etomica.AtomSequencerSimple;
import etomica.AtomsetIterator;
import etomica.IteratorFactory;
import etomica.Species;
import etomica.AtomSequencer.Factory;

public class IteratorFactoryCell extends IteratorFactory {
    
    public static final IteratorFactoryCell INSTANCE = new IteratorFactoryCell();
   
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
        ApiInnerVariable pairIterator = new ApiInnerVariable(
                    new AtomIteratorBasis(),
                    new AtomIteratorNbrCell(species[0],true));
        return new ApiIntragroup(pairIterator);
    }
    public Factory moleculeSequencerFactory() {
        return AtomSequencerSimple.FACTORY;
    }
}