package simulate;

public class SpeciesWalls extends Species {

    int borderTol;
    boolean boundary;

    public SpeciesWalls() {
        this(1,1);
    }
    
    public SpeciesWalls(int n, int na) {
        super(n, na, new AtomWall());
        colorScheme.setBaseColor(Constants.DARK_RED);
        add(new ConfigurationMoleculeWallsParallel());
    }

  public void initializeMolecules() {
    setThickness(4);
  }

  public Molecule makeMolecule() {
    return new MoleculeWall(this,nAtomsPerMolecule);
  }

  public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
//  public void initializeSpecies(Phase phase) {
//    configurationMolecule.initializeCoordinates();
//  }
}
