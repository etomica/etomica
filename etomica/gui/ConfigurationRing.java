package simulate;
import java.util.*;
import java.awt.*;

/**
 * Fills phase with molecules on a lattice, taking each molecule in successive order
 * from the linked list of molecules.  Takes no special action when list moves from
 * one species to the next.
 * Wall "molecules" are ignored becuase the (super) add method will not add them.
 *
 * Need to improve this to handle different dimensions more elegantly
 */
 
public class ConfigurationRing extends Configuration {

    private boolean fill;
    
    public ConfigurationRing() {
        super();
        setFillVertical(true);
    }
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    public void initializeCoordinates() {
        if(parentPhase == null) {return;}
        
        double Lx = parentPhase.dimensions().component(0);
        double Ly = 0.0;
        if(Simulation.space().D()>1)  Ly = parentPhase.dimensions().component(1);

    // Count number of molecules
        int sumOfMolecules = 0;
        for(int j=0; j<species.size(); j++) { 
            Species.Agent s = (Species.Agent)species.elementAt(j);
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            sumOfMolecules += s.moleculeCount();
        }
        
        if(sumOfMolecules == 0) {return;}
        
        Space.Vector[] rLat;
        if(Simulation.space.D() == 1) {
            rLat = lineLattice(sumOfMolecules, Lx);
        }
        else {
            rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
        }
        
   // Place molecules     
        int i = 0;
        for(int j=0; j<species.size(); j++) {
            Species.Agent s = (Species.Agent)species.elementAt(j);
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            for(Molecule m=s.firstMolecule(); m!=s.terminationMolecule(); m=m.nextMolecule()) {
                m.setCOM(rLat[i]);
//                m.firstAtom().r.setRandom(1.0);
                i++;
            }
        }
//        parentPhase.parentSimulation.space().clearCells();
 //       for(Atom a=parentPhase.firstAtom(); a!=null; a=a.nextAtom()) {
 //           a.coordinate.assignCell();}
        initializeMomenta();
    }
}
