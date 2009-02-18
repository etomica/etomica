package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;

public class BoxMoleculeEvent extends BoxEvent implements IBoxMoleculeEvent {
        
        public BoxMoleculeEvent(IBox box, IMolecule mole) {
            super(box);
            this.molecule = mole;
        }

        /* (non-Javadoc)
         * @see etomica.box.IBoxAtomEvent#getAtom()
         */
        public IMolecule getMolecule() {
            return molecule;
        }
        
        private final IMolecule molecule;
        private static final long serialVersionUID = 1L;
}
