package etomica.graphics;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    
    public ColorSchemeRandom(Simulation sim) {
        super(sim);
    }
    
    public void colorAllAtoms(Phase phase) {
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if(atomColors[a.getGlobalIndex()] == null) {
                atomColors[a.getGlobalIndex()] = ConstantsGraphic.randomColor();
            }
        }
    }
}