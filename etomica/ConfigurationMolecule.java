package simulate;
import java.util.*;
import java.awt.*;

/**
 * General class for assignment of coordinates to all atoms
 * Places atoms of each molecule in species in same position with respect to the molecule's center of mass
 */

public abstract class ConfigurationMolecule extends Component {
  
  Species parentSpecies;
  protected final double[] dim = new double[Space.D];
  
  public ConfigurationMolecule(){
    }
    
/*  public void addNotify() {
    super.addNotify();
    if(getParent() instanceof Species) {
       parentSpecies = (Species)getParent();
    }
    else {  //use an exception here?
        System.out.println("Error:  ConfigurationMolecule should be added only to a Species");
    }
  }
*/  
  public void initializeCoordinates() {
    if(parentSpecies == null) {return;}
    for(Molecule m=parentSpecies.firstMolecule(); m!=parentSpecies.terminationMolecule(); m=m.getNextMolecule()) {
        initializeCoordinates(m);
    }
    computeDimensions();
  }
  
  public double[] moleculeDimensions() {return dim;}
  
  public abstract void initializeCoordinates(Molecule m);
  protected abstract void computeDimensions();
  
}
