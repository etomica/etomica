package etomica.nbr.cell;

import etomica.atom.AtomSequencerFactory;
import etomica.atom.iterator.IteratorFactorySimple;

//TODO revise after consideration of PotentialMaster

public class IteratorFactoryCell extends IteratorFactorySimple {

    public AtomSequencerFactory interactionAtomSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomSequencerFactory interactionMoleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
    public AtomSequencerFactory moleculeSequencerFactory() {
        return AtomSequencerCell.FACTORY;
    }
}