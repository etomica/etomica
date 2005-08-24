package etomica.graphics;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.phase.Phase;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorListTabbed iterator = new AtomIteratorListTabbed();
    
    public ColorSchemeRandom() {
        super();
    }
    
    public void colorAllAtoms(Phase phase) {
        iterator.setList(phase.getSpeciesMaster().atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if(a.allatomAgents[agentIndex] == null) {
                a.allatomAgents[agentIndex] = ConstantsGraphic.randomColor();
            }
        }
    }
}//end of ColorSchemeRandom