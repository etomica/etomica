package simulate;
import java.util.*;
import java.awt.*;

/**
 * General class for assignment of coordinates to all atoms
 * Places atoms of each molecule in species in same position with respect to the molecule's center of mass
 */

public abstract class ConfigurationMolecule extends Component {
  
  protected Species parentSpecies;  //some subclasses may want to take an action on setting species, so don't make public
  protected final double[] dim = new double[Simulation.D];
  
  public ConfigurationMolecule(){
    }
  
  public void initializeCoordinates() {
    if(parentSpecies == null) {return;}
    Enumeration e = parentSpecies.agents.elements();
    while(e.hasMoreElements()) {
        Species.Agent agent = (Species.Agent)e.nextElement();
        for(Molecule m=agent.firstMolecule(); m!=agent.terminationMolecule(); m=m.nextMolecule()) {
            initializeCoordinates(m);
        }
    }
    computeDimensions();
  }
        
        
  public void initializeCoordinates(Phase p) {
    if(parentSpecies == null) {return;}
    Species.Agent agent = parentSpecies.getAgent(p);
    for(Molecule m=agent.firstMolecule(); m!=agent.terminationMolecule(); m=m.nextMolecule()) {
        initializeCoordinates(m);
    }
    computeDimensions();
  }
  
  public void setParentSpecies(Species s) {parentSpecies = s;}
  public Species parentSpecies() {return parentSpecies;}
  
  public double[] moleculeDimensions() {return dim;}
  
  public abstract void initializeCoordinates(Molecule m);
  protected abstract void computeDimensions();
  
}
