package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDiskWell extends SpeciesDisks {
  double lambda;
  Color wellColor;

  public SpeciesDiskWell(Simulation ps) {
    this(ps, 20,1);
  }
  public SpeciesDiskWell(Simulation ps, int nM, int nA) {
    super(ps, nM, nA, new AtomType.Well(1.0,Color.black,0.1,1.5));
  }
  
  public SpeciesDiskWell() {
        setSpeciesIndex(0);
        setNAtomsPerMolecule(1);
        protoType = new AtomType.Well(1.0,Color.black,0.1,1.5);  //mass, color, diameter, lambda
        colorScheme = new ColorSchemeByType();
        configurationMolecule = new ConfigurationMoleculeLinear();
        setNMolecules(20);
  }
  
  public double getLambda() {return ((AtomType.Well)protoType).lambda();}
  public void setLambda(double lam) {((AtomType.Well)protoType).setLambda(lam);}
}