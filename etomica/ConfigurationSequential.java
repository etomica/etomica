package etomica;
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
    
    public ConfigurationSequential(Space space) {
        super(space);
        setFillVertical(true);
    }
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    //need to revise to use given argument instead of parentphase
    public void initializeCoordinates(Atom atom) {
        AtomGroup atoms;
        
        if(!(atom instanceof AtomGroup)) {
            atom.coord.position().E(0.0);
            return;
        }
        else atoms = (AtomGroup)atom;
            
        Phase phase = atoms.parentPhase();
        
        
        if(phase == null) return;
        double Lx = phase.dimensions().component(0);
        double Ly = 0.0;
        if(phase.parentSimulation().space().D()>1)  Ly = phase.dimensions().component(1);

    // Count number of molecules
        int sumOfMolecules = 0;
        for(SpeciesAgent s=phase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            sumOfMolecules += s.moleculeCount();
        }
        
        if(sumOfMolecules == 0) {return;}
        
        Space.Vector[] rLat;
        if(phase.parentSimulation().space().D() == 1) {
            rLat = null; //for redesign
// commented for redesign            rLat = lineLattice(sumOfMolecules, Lx);
        }
        else {
            rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
        }
        
   // Place molecules     
        int i = 0;
        for(SpeciesAgent s=phase.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.parentSpecies() instanceof SpeciesWalls) {continue;}
            AtomIterator iterator = s.childIterator;
            iterator.reset();
            while(iterator.hasNext()) {
                iterator.next().coord.translateTo(rLat[i]);
                i++;
            }
        }
   //     initializeMomenta(phase.speciesMaster());
        initializeMomenta(atom);
    }
}
