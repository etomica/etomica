package simulate;

import java.awt.Graphics;

/** Four walls arranged as a piston-cylinder apparatus.  All but one
 *  wall (the piston) is stationary.
 */
 
public class SpeciesPistonCylinder extends Species {

    private int orientation;
    private int thickness;
    private transient final int[] sideWallOrigin = new int[Space.D];
    
    public SpeciesPistonCylinder() {
        this(1,4);  //1 molecule, 4 atoms
    }
    
    public SpeciesPistonCylinder(int n, int na) {
        super(1, 4, new AtomHardWall(null));
        colorScheme.setBaseColor(Constants.DARK_RED);
        this.add(new ConfigurationMoleculeWallsBoundary());
        setOrientation(Constants.NORTH);
        initializeMolecules();
    }
    
    public void setNMolecules(int i) {super.setNMolecules(1);}  //override so nMolecules cannot be changed
    public void setNAtomsPerMolecule(int i) {super.setNAtomsPerMolecule(4);} //likewise

  public void initializeMolecules() {
    setThickness(4);
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {a.setStationary(true);}
    firstAtom().setStationary(false);  //first atom is the piston
    ((AtomWall)firstAtom().getNextAtom()).setLength(0.1*Double.MAX_VALUE);
    ((AtomWall)lastAtom()).setLength(0.1*Double.MAX_VALUE);
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
    
  public final int getThickness() {return thickness;}
  public final void setThickness(int t) {
    thickness = t;
    Atom terminator = terminationAtom();
    for(Atom a=firstAtom(); a!=terminator; a=a.getNextAtom()) {
        ((AtomWall)a).setThickness(t);
    }
  } 
  
  //override to make sure full side walls are drawn when scaling image
  public void draw(Graphics g, int[] origin, double scale) {

    double toPixels = scale*Phase.TO_PIXELS;
    Atom nextSpeciesAtom = lastAtom().getNextAtom();
    Molecule last = lastMolecule.getNextMolecule();
    if(scale != 1.0) {
        int newThickness = Math.max(1,(int)(scale*thickness));
        for(Atom a=firstAtom(); a!=nextSpeciesAtom; a=a.getNextAtom()) {
          ((AtomWall)a).setThickness(newThickness);
        }
    }
    
    Atom a;
    
    //draw piston
    a = firstAtom();
    colorScheme.setAtomColor(a);
    a.draw(g,origin,scale);
    
    //draw bottom
    a = firstAtom().getNextAtom().getNextAtom();
    colorScheme.setAtomColor(a);
    a.draw(g,origin,scale);

    //draw side walls
    if(orientation == Constants.NORTH || orientation == Constants.SOUTH) {
        sideWallOrigin[0] = origin[0];
        sideWallOrigin[1] = 0;
    }
    else {
        sideWallOrigin[0] = 0;
        sideWallOrigin[1] = origin[1];
    }
    a = firstAtom().getNextAtom();
    colorScheme.setAtomColor(a);
    a.draw(g,sideWallOrigin,scale);

    a = lastAtom();
    colorScheme.setAtomColor(a);
    a.draw(g,sideWallOrigin,scale);

    //reset thickness to original value
    setThickness(thickness);
  }
}
