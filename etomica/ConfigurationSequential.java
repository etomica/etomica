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
    private final AtomIteratorSequential iterator = new AtomIteratorSequential();
    
    public ConfigurationSequential(Space space) {
        super(space);
        setFillVertical(true);
    }
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    public void initializeCoordinates(AtomIterator iterator) {
        
        Phase phase = iterator.getBasis().parentPhase();
        
        if(phase == null) return;
        double Lx = phase.dimensions().component(0);
        double Ly = 0.0;
        double Lz = 0.0;
        if(phase.parentSimulation().space().D()>1)  Ly = phase.dimensions().component(1);
        if(phase.parentSimulation().space().D()>2)  Lz = phase.dimensions().component(2);

        int sumOfMolecules = iterator.size();
        
        if(sumOfMolecules == 0) return;
        
        Space.Vector[] rLat;
        switch(space.D()) {
            case 1:
                rLat = lineLattice(sumOfMolecules, Lx);
                break;
            default:
            case 2:
                rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
                break;
            case 3:
                rLat = ConfigurationFcc.lattice(sumOfMolecules);
                break;
        }
        
   // Place molecules     
        int i = 0;
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.next();
            if(a.parentSpecies() instanceof SpeciesWalls) continue;
            a.coord.translateTo(rLat[i]);
            i++;
        }
   //     initializeMomenta(phase.speciesMaster());
        initializeMomenta(iterator.getBasis());
    }
}
