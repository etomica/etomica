package etomica.nbr.cell;

import etomica.ApiMolecule;
import etomica.AtomSequencerFactory;
import etomica.AtomsetIteratorMolecule;
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
    public AtomsetIteratorMolecule makeInterspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiInterspecies1ACell(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiInterspeciesAACell(D, species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiIntraspecies1ACell(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAACell(D, species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}