/*
 * Created on Nov 15, 2004
 *
 */
package etomica.models.hexane;

import etomica.AtomFactory;
import etomica.AtomIndexManager;
import etomica.Space;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomSequencerFactory;

/**
 * @author nancycribbin
 *
 * Factory which creates hexane molecules based on Dr. Monson's data.  The molecule has
 * 6 carbon atoms as child atoms, but ignores the existance of hydrogen atoms.
 * 
 */
public class AtomFactoryHexane extends AtomFactoryHomo {
    
    public AtomFactoryHexane(Space space, AtomSequencerFactory sequencerFactory,
            AtomIndexManager indexManager, AtomFactory factory){
        super(space, sequencerFactory, indexManager, factory, 6,
            new ConformationHexane(space));
    }
}