package simulate;

public class SpeciesWalls extends Species {

    public static final int NORTH = 0;
    public static final int SOUTH = 1;
    public static final int EAST = 2;
    public static final int WEST = 3;
    
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

  public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
}
