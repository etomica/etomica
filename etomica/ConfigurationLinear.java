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
    public void initializeCoordinates(Atom group) {
        
        if(!(group instanceof AtomGroup) 
            || ((AtomGroup)group).childCount() < 2) return; //nothing to do if group is just one atom
            
        work.E(group.coord.position());//save original COM position
        double xNext = 0.0;
        for(Atom a=((AtomGroup)group).firstChild(); a!=null; a=a.nextAtom()) {
            a.coord.translateTo(work);  //put all atoms at same point
            a.coord.translateBy(xNext,orientation);  //move xNext distance in direction orientation
            xNext += bondLength;
        }
        group.coord.translateTo(work);  //shift molecule to original COM
    }//end of initializeCoordinates
}//end of ConfigurationLinear
      
