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
 
public class ConfigurationSequential extends Configuration {

    private boolean fill;
    
    public ConfigurationSequential() {
        super();
        setFillVertical(true);
    }
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    public void initializeCoordinates(Phase phase) {
        
        double Lx = phase.dimensions().component(0);
        double Ly = 0.0;
        if(phase.parentSimulation().space().D()>1)  Ly = phase.dimensions().component(1);

    // Count number of molecules
        int sumOfMolecules = 0;
        for(Species.Agent s=phase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            sumOfMolecules += s.moleculeCount();
        }
        
        if(sumOfMolecules == 0) {return;}
        
        Space.Vector[] rLat;
        if(phase.parentSimulation().space().D() == 1) {
            rLat = lineLattice(sumOfMolecules, Lx);
        }
        else {
            rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
        }
        
   // Place molecules     
        int i = 0;
        for(Species.Agent s=phase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            for(Molecule m=s.firstMolecule(); m!=s.terminationMolecule(); m=m.nextMolecule()) {
                m.setCOM(rLat[i]);
                i++;
            }
        }
        initializeMomenta(phase);
    }
}
