package simulate;
import java.util.*;
import java.awt.*;

/**
 * Fills phase with molecules on a lattice, taking each molecule in successive order
 * from the linked list of molecules.  Takes no special action when list moves from
 * one species to the next.
 * Wall "molecules" are ignored becuase the (super) add method will not add them.
 */
 
public class ConfigurationSequential extends Configuration {

    private boolean fill;
    
    public ConfigurationSequential() {
        super();
        setFillVertical(true);
    }
    
    public void setFillVertical(boolean b) {fill = b;}
    public boolean getFillVertical() {return fill;}
    
    public void initializeCoordinates() {
        if(parentPhase == null) {return;}
        
        double Lx = parentPhase.getBounds().width/Phase.TO_PIXELS;
        double Ly = parentPhase.getBounds().height/Phase.TO_PIXELS;

        int sumOfMolecules = 0;
        for(int j=0; j<species.size(); j++) {   
            Species s = (Species)species.elementAt(j);
            sumOfMolecules += s.getNMolecules();
        }
        
        double[][] rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
        
        int i = 0;
        for(int j=0; j<species.size(); j++) {
            Species s = (Species)species.elementAt(j);
            for(Molecule m=s.firstMolecule(); m!=s.terminationMolecule(); m=m.getNextMolecule()) {
                m.setCOM(rLat[i]);
                i++;
            }
        }
        initializeMomenta();
    }
}
