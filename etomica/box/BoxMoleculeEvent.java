package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;

public class BoxMoleculeEvent extends BoxEvent implements IBoxMoleculeEvent {
        
        public BoxMoleculeEvent(IBox box, IMolecule mole) {
            super(box);
            this.molecule = mole;
        }
        
        public IMolecule getMolecule() {
            return molecule;
        }

        
        protected IMolecule molecule = null;
        private static final long serialVersionUID = 1L;
}
