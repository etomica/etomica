package etomica.graphics;
import etomica.*;
import etomica.atom.iterator.AtomIteratorListSimple;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorListSimple iterator = new AtomIteratorListSimple();
    
    public ColorSchemeRandom() {
        super();
    }
    
    public void colorAllAtoms(Phase phase) {
        iterator.setList(phase.speciesMaster.atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if(a.allatomAgents[agentIndex] == null) {
                a.allatomAgents[agentIndex] = ConstantsGraphic.randomColor();
            }
        }
    }
}//end of ColorSchemeRandom