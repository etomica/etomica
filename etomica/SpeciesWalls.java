package simulate;

public class SpeciesWalls extends Species {

    int borderTol;
    boolean boundary;

    public SpeciesWalls() {
        this(1,1);
    }
    
    public SpeciesWalls(int n, int na) {
        super(n, na, new AtomHardWall(null));
        colorScheme.setBaseColor(Constants.DARK_RED);
        this.add(new ConfigurationMoleculeWallsParallel());
    }

  public void initializeMolecules() {
    setThickness(4);
  }

  public final int getThickness() {return ((AtomHardWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomHardWall)firstAtom()).setThickness(t);}
    
//  public void initializeSpecies(Phase phase) {
//    configurationMolecule.initializeCoordinates();
//  }
}
