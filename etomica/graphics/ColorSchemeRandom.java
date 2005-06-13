package etomica.graphics;
import etomica.Atom;
import etomica.Phase;
import etomica.atom.iterator.AtomIteratorListTabbed;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorListTabbed iterator = new AtomIteratorListTabbed();
    
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