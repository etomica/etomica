package simulate;

/** Four walls arranged as a piston-cylinder apparatus.  All but one
 *  wall (the piston) is stationary.
 */
 
public class SpeciesPistonCylinder extends Species {

    private int orientation;
    
    public SpeciesPistonCylinder() {
        this(1,4);  //1 molecule, 4 atoms
    }
    
    public SpeciesPistonCylinder(int n, int na) {
        super(1, 4, new AtomHardWall(null));
        colorScheme.setBaseColor(Constants.DARK_RED);
        this.add(new ConfigurationMoleculeWallsBoundary());
        setOrientation(Constants.NORTH);
        for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {a.setStationary(true);}
        firstAtom().setStationary(false);  //first atom is the piston
    }
    
    public void setNMolecules(int i) {super.setNMolecules(1);}  //override so nMolecules cannot be changed
    public void setNAtomsPerMolecule(int i) {super.setNAtomsPerMolecule(4);} //likewise

  public void initializeMolecules() {
    setThickness(4);
  }
  
  public int getOrientation() {return orientation;}
  public void setOrientation(int i) {
    orientation = i;
    ((ConfigurationMoleculeWallsBoundary)configurationMolecule).setFirstWallPosition(orientation);
    switch(orientation) {
        case Constants.NORTH:
        case Constants.SOUTH:
        case Constants.EAST:
        case Constants.WEST:
    }
  }
    
  public final int getThickness() {return ((AtomHardWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomHardWall)firstAtom()).setThickness(t);}
    
}
