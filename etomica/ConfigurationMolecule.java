package simulate;
import java.util.*;
import java.awt.*;

/**
 * General class for assignment of coordinates to all atoms
 * Places atoms of each molecule in species in same position with respect to the molecule's center of mass
 */

public abstract class ConfigurationMolecule extends Component {
  
  Species parentSpecies;
  protected final double[] dim = new double[Simulation.D];
  
  public ConfigurationMolecule(){
    }
    
  public void initializeCoordinates() {
    if(parentSpecies == null) {return;}
    for(Molecule m=parentSpecies.firstMolecule(); m!=parentSpecies.terminationMolecule(); m=m.nextMolecule()) {
        initializeCoordinates(m);
    }
    computeDimensions();
  }
  
  public double[] moleculeDimensions() {return dim;}
  
  public abstract void initializeCoordinates(Molecule m);
  protected abstract void computeDimensions();
  
}
