package etomica.graphics;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class ColorSchemeRandom extends ColorSchemeCollective {
    
    private final AtomIteratorListTabbed iterator = new AtomIteratorListTabbed();
    
    public ColorSchemeRandom(Simulation sim) {
        super(sim);
    }
    
    public void colorAllAtoms(Phase phase) {
        iterator.setList(phase.getSpeciesMaster().atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if(atomColors[a.getGlobalIndex()] == null) {
                atomColors[a.getGlobalIndex()] = ConstantsGraphic.randomColor();
            }
        }
    }
}