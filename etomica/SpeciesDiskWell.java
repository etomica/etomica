package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDiskWell extends SpeciesDisks {
  double lambda;
  Color wellColor;

  public SpeciesDiskWell(int nM, int nA) {
    this(Simulation.instance, nM, nA);
  }
  public SpeciesDiskWell(Simulation sim, int nM, int nA) {
    super(sim, nM, nA, new AtomType.Well(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE,1.5));
  }
  
  public SpeciesDiskWell() {
    this(Simulation.instance);
  }
    
  public SpeciesDiskWell(Simulation sim) {
        super(sim);
        setAtomsPerMolecule(1);
        protoType = new AtomType.Well(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE,1.5);  //mass, color, diameter, lambda
        moleculeConfiguration = new Molecule.Configuration.Linear(this);
        setNMolecules(Default.MOLECULE_COUNT);
  }
  
  public double getLambda() {return ((AtomType.Well)protoType).lambda();}
  public void setLambda(double lam) {((AtomType.Well)protoType).setLambda(lam);}
}