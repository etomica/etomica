package etomica;
import etomica.units.Dimension;

/**
 * Places atoms in a straight line. 
 *
 * @author David Kofke
 */

public class ConfigurationLinear extends Configuration {
          
    private double bondLength = 0.5*Default.ATOM_SIZE;
    private Space.Vector orientation;
    private double[] angle;
    
    public ConfigurationLinear(Space space) {
        super(space);
        orientation = space.makeVector();
        angle = new double[space.D()];
        setAngle(0,etomica.units.Degree.UNIT.toSim(45.));
    }
              
    public void setBondLength(double b) {
        bondLength = b;
    }
    public double getBondLength() {return bondLength;}
    public Dimension getBondLengthDimension() {return Dimension.LENGTH;}
    
    //need to re-express this in terms of a Space.Orientation object
    public void setAngle(int i, double t) {//t in radians
        angle[i] = t;
        switch(angle.length) {
            case 1:
                return;
            case 2:
                setOrientation(new Space2D.Vector(Math.cos(angle[0]),Math.sin(angle[0])));
                return;
//            case 3:
//                setOrientation(new Space3D.Vector(Math.sin(angle[1])*Math.cos(angle[0]),
//                                                  Math.sin(angle[1])*Math.sin(angle[0]),
//                                                  Math.cos(angle[1])));
//                return;
        }
    }
    public double getAngle(int i) {return angle[i];}
    public void setOrientation(Space.Vector e) {orientation.E(e);}
              
    /**
    * Sets all atoms coordinates to lie on a straight line along the x-axis, with the
    * center of mass unchanged from the value before method was called
    */
    public void initializeCoordinates(AtomIterator iterator) {
        
        int size = iterator.size();
        if(iterator.size() == 0) return;
            
        double xNext = -bondLength*0.5*(double)(size-1);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.next();
            a.coord.translateTo(space.origin());
            a.coord.translateBy(xNext,orientation);  //move xNext distance in direction orientation
            xNext += bondLength;
        }
    }//end of initializeCoordinates
}//end of ConfigurationLinear
      
