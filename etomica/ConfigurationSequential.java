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
    private Space.Vector dimensions;
    
    public ConfigurationSequential(Space space) {
        super(space);
        setFillVertical(true);
        dimensions = space.makeVector();
        dimensions.E(Default.BOX_SIZE);
    }
    
    public void setDimensions(Space.Vector dimensions) {this.dimensions = dimensions;}
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    public void initializePositions(AtomIterator[] iterators) {

        AtomIteratorCompound iterator = new AtomIteratorCompound(iterators);//lump 'em all together

        double Lx = dimensions.component(0);
        double Ly = 0.0;
        double Lz = 0.0;
        if(dimensions.length()>1)  Ly = dimensions.component(1);
        if(dimensions.length()>2)  Lz = dimensions.component(2);

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
    }
}
