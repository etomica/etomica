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
    public AtomsetIteratorMolecule makeInterSpeciesPairIterator(Species[] species) {
        //TODO make the 1A iterators
        AtomsetIteratorMolecule api1A = new ApiInterspecies1ACell(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiInterspeciesAACell(D, species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraSpeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = null;//new ApiIntraSpecies1ACell(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAACell(D, species);
        return new ApiMolecule(api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}