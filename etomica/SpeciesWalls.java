package simulate;
import java.awt.*;//for Graphics
import java.beans.Beans;

public class SpeciesWalls extends Species {

    int borderTol;
    boolean boundary;

    public SpeciesWalls() {
        this(1,1);
    }
    
    public SpeciesWalls(int n, int na) {
        super(n,na);
        configurationMolecule = new ConfigurationMoleculeWallsParallel();
        initializeMolecules();
    }

  public void initializeMolecules() {
    setThickness(4);
    colorScheme.setBaseColor(Constants.DARK_RED);
  }

  public Molecule makeMolecule() {
    return new MoleculeWall(this,nAtomsPerMolecule);
  }

  public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
  public void addNotify() {
      super.addNotify();
  }
    
  public void initializeSpecies(Phase phase) {
    configurationMolecule.initialize();
  }
}
