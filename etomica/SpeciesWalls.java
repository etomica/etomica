package simulate;
import java.awt.Color;

public class SpeciesWalls extends Species {

    public static final int NORTH = 0;
    public static final int SOUTH = 1;
    public static final int EAST = 2;
    public static final int WEST = 3;

/** 
 *  Wall type array.  Each atom has its own type, which specifies its length and orientation.
 *  Examples of use:
 *  Use one molecule, with atoms of different type to set up box, cylinder, etc.
 *  Use several molecules each with one atom to set up parallel walls.
 */
    public AtomType.Wall[] protoType;

    /**
    * Default constructor.  Creates species containing 1 molecule with 1 horizontal wall atom.
    */
    public SpeciesWalls() {
        this(1,1,Double.MAX_VALUE,0);
    }
    
    public SpeciesWalls(int nM, int nA, double length, int angle) {  //angle is in degrees
        protoType = new AtomType.Wall[nA];
        for(int i=0; i<nA; i++) {protoType[i] = new AtomType.Wall(1.0, Color.black, length, angle);}  // arguments are mass, color, length, angle(degrees)
 //       this(nM, protoType);  //can't do this because must be first line in constructor
        setSpeciesIndex(1);     
        atomsPerMolecule = nA;
        setNMolecules(nM);
        
        colorScheme = new ColorSchemeNull();
        this.add(new ConfigurationMoleculeWallsParallel());
    }

    public SpeciesWalls(int nM, AtomType.Wall[] type) {  
        setSpeciesIndex(1);
        protoType = type;
        atomsPerMolecule = type.length;
 //       setNMolecules(nM);
        
        colorScheme = new ColorSchemeNull();
        this.add(new ConfigurationMoleculeWallsParallel());
    }

    protected Molecule makeMolecule(Phase phase) {
        return new Molecule(this, phase, protoType);
    } 
              
    // Exposed Properties --- not implemented now because they must tie to AtomType array
    
/*  public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
  public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
                
    public final double getLength() {return protoType.length();}
    public void setLength(double d) {protoType.setLength(d);}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}*/
}
