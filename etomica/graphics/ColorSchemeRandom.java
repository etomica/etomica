package etomica.graphics;
import etomica.*;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorListSimple iterator = new AtomIteratorListSimple();
    
    public ColorSchemeRandom() {
        super();
    }
    
    public void colorAllAtoms(Phase phase) {
        iterator.setList(phase.speciesMaster.atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.next();
            if(a.allatomAgents[agentIndex] == null) {
                a.allatomAgents[agentIndex] = ConstantsGraphic.randomColor();
            }
        }
    }
}//end of ColorSchemeRandom