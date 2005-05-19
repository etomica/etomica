package etomica.nbr.cell;

import etomica.IteratorFactory;
import etomica.Species;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.iterator.ApiMolecule;
import etomica.atom.iterator.ApiSpecies11;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;

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
        ApiSpecies11 api11 = new ApiSpecies11(species);
        return new ApiMolecule(api11, api1A, apiAA);
    }
    public AtomsetIteratorMolecule makeIntraspeciesPairIterator(Species[] species) {
        AtomsetIteratorMolecule api1A = new ApiIntraspecies1ACell(D, species);
        AtomsetIteratorPhaseDependent apiAA = new ApiIntraspeciesAACell(D, species);
        ApiSpecies11 api11 = new ApiSpecies11(species);
        return new ApiMolecule(api11, api1A, apiAA);
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}