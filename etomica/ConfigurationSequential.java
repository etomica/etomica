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
        if(parentPhaseSpace == null) {return;}
        
        double Lx = parentPhaseSpace.getBounds().width/DisplayConfiguration.SIM2PIXELS;
        double Ly = parentPhaseSpace.getBounds().height/DisplayConfiguration.SIM2PIXELS;

        int sumOfMolecules = 0;
        for(int j=0; j<species.size(); j++) {   
            Species s = (Species)species.elementAt(j);
            sumOfMolecules += s.getNMolecules();
        }
        
        PhaseSpace2D.Vector[]  rLat = squareLattice(sumOfMolecules, Lx, Ly, fill); 
        
        int i = 0;
        for(int j=0; j<species.size(); j++) {
            Species s = (Species)species.elementAt(j);
            for(Molecule m=s.firstMolecule(); m!=s.terminationMolecule(); m=m.nextMolecule()) {
                m.setCOM(rLat[i]);
                i++;
            }
        }
        initializeMomenta();
    }
}
