package etomica;

/* History of changes
 * 07/29/02 (DAK) In SequentialIterator, changed reset(IteratorDirective) to handle null atom1()
 * 08/13/02 (DAK) In SequentialIterator, changed to allow basis to be leaf atom.  This introduces
 *               a singlet iterator and currentIterator to handle the general case.
 * 12/05/02 (DAK) Revised reset(IteratorDirective) so that it properly handles a reference atom
 *                this is descended from the basis, but not is a direct child of the basis.
 */

public class IteratorFactorySimple implements IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
        
    public AtomSequencer makeSimpleSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
        
    public AtomsetIterator makeIntergroupNbrPairIterator() {return new ApiIntergroup();}
    public AtomsetIterator makeIntragroupNbrPairIterator() {return new ApiIntragroup();}
    
    public AtomSequencer makeNeighborSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
    
    public Class simpleSequencerClass() {return AtomSequencerSimple.class;}
    
    public Class neighborSequencerClass() {return AtomSequencerSimple.class;}
    
    public AtomSequencer.Factory simpleSequencerFactory() {return AtomSequencerSimple.FACTORY;}
    
    public AtomSequencer.Factory neighborSequencerFactory() {return AtomSequencerSimple.FACTORY;}
   
}