package etomica;
import etomica.units.*;
import etomica.electrostatics.*;

/**
 * AtomType holds atom parameters.  
 * It is used to set the general features of the atom (e.g., whether it is a wall, sphere, etc.),
 * and to prescribe how it is drawn to the screen.
 * AtomType is responsible for selecting the appropriate type of coordinate needed to describe the
 * position and (perhaps) orientation of the atom.
 * The AtomType of an atom is set by Species when it constructs a molecule.  Each Atom has an instance variable
 * named "type" that holds the AtomType object; this may be different for each atom in a molecule, or
 * it may refer to a common AtomType object, as prescribed by the Species.
 * AtomType could also be used to define particular elemental atoms (Carbon, Oxygen, etc.).
 * 
 */

public abstract class AtomType implements java.io.Serializable {
    public static String getVersion() {return "AtomType:01.11.20";}
    private double mass, rm;
    private ElectroType electroType;
    private String name;
    private /*final*/ AtomFactory creator;
    
    public AtomType(AtomFactory creator) {
        this(creator, Default.ATOM_MASS);
    }
    
    public AtomType(AtomFactory creator, double m) {
        this.creator = creator;
        setMass(m);
        setElectroType(new ElectroType.Null());
    }
    
    public AtomFactory creator() {return creator;}
    /**
     * Returns default coordinate type, which has no orientational component.
     * Override for atom types that require other coordinate features.
     */
/*    public Space.Coordinate makeCoordinate(Atom a) {
        return a.parentSimulation().space().makeCoordinate(a);
    }
*/
    public void initialize(Atom a) {
        a.coord.setMass(mass);
    }
        
    public final void setElectroType(ElectroType et) {
        electroType = et;
    }
    public final ElectroType electroType() {return electroType;}
            
    /**
    * Sets  mass of this atom and updates reciprocal mass accordingly.  Setting
    * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass 
    * to be set to zero.
    * 
    * @param mass   new value for mass
    */
    public void setMass(double m) {
        mass = m;
        rm = (mass==Double.MAX_VALUE) ? 0.0 : 1.0/mass;
    }
    public final double getMass() {return mass;}
    public final double mass() {return mass;}
    public final Dimension getMassDimension() {return Dimension.MASS;}
    public final void setRm(double rm) {
        setMass( (rm==0) ? Double.MAX_VALUE : 1.0/rm);
    }
    public final double rm() {return rm;}


    // Sphere-shaped atom.  Assumed to be in 2D
    public static class Sphere extends AtomType {
        
        double diameter, radius;
        
        public Sphere(AtomFactory creator) {
            this(creator, Default.ATOM_MASS, Default.ATOM_SIZE);
        }
        public Sphere(AtomFactory creator, double m, double d) {
            super(creator, m);
            setDiameter(d);
        }
                    
        public final double diameter() {return diameter;}
        public final double radius() {return radius;}
        
        /**
        * Sets diameter of this atom and updates radius accordingly.
        *
        * @param d   new value for diameter
        */
        public void setDiameter(double d) {diameter = d; radius = 0.5*d;}
    }
    
    /**
     * Atom type for a sphere that has some feature depending upon an orientation coordinate.
     * For example an orientational dependent potential may be attached to an otherwise spherical atom
     */
    public final static class OrientedSphere extends Sphere implements SphericalTop {
        
        private final double[] I = new double[3];
        public OrientedSphere(AtomFactory creator, double m, double d) {
            super(creator,m,d);
            updateI();
        }
        public double[] momentOfInertia() {return I;}
        
        public Space.Coordinate makeCoordinate(Atom a) {
            return a.parentSimulation().space().makeCoordinate(a); //override changes nothing, but this may change if revise method in Space
        }
        
        private void updateI() {
            if(I == null) return;
            I[0] = 0.4*this.mass()*radius()*radius();  //moment of inertia of a sphere = 2/5 m R^2 (should modify to arbitrary dimension)
            I[1] = I[0];
            I[2] = I[1];
        }
        
        public void setMass(double m) {
            super.setMass(m);
            updateI();
        }
        public void setDiameter(double d) {
            super.setDiameter(d);
            updateI();
        }
    }
    
    // Sphere with a concentric well.
    public final static class Well extends Sphere {  
        
        private double lambda;                    //diameter of well, in units of core diameter
        private double wellDiameter, wellRadius;  //size of well, in simulation units
        
        public Well(AtomFactory creator, double m, double d, double l) {
            super(creator, m, d);
            setDiameter(d);
            setLambda(l);
        }
                    
        public final double lambda() {return lambda;}
        public final double wellDiameter() {return wellDiameter;}
        public final double wellRadius() {return wellRadius;}
        
        public final void setDiameter(double d) {
            super.setDiameter(d); 
            setLambda(lambda);
        }
        public final void setLambda(double l) {
            lambda = l; 
            wellDiameter = lambda*diameter(); 
            wellRadius = 0.5*wellDiameter;
        }
    }

    // Wall-shaped atom.  Arbitrary dimension.
    public final static class Wall extends AtomType {
        
        int thickness = 4;  //thickness when drawn to screen (if horizontal or vertical)
        private boolean vertical, horizontal, wide;
        private double cosX, sinX, tanX, cosY, sinY, tanY, cosZ, sinZ, tanZ;
//        double[] f = new double[Space.D];   //force on wall
//        double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
        protected double xAngle, yAngle, zAngle;   //orientation (radians) with respect to positive x-axis, taking r at the origin
        protected double length;  //length of line in simulation units
        protected boolean longWall;  //set to true if wall is infinite in length
        //  these may belong in integratorAgent
        protected double temperature = Kelvin.UNIT.toSim(300.);
        protected boolean adiabatic = true;
        private Constants.Alignment alignment;
        public double pAccumulator, qAccumulator; //net sum of momentum and heat transferred to wall
        
        public Wall(AtomFactory creator, double m, double l, double x, double y, double z) {
            super(creator, m);
            setLength(l);
            setXAngle(x);
            setYAngle(y);
            setZAngle(z);
        }
        
        private void checkAlignment() {
          horizontal = (xAngle == 0.0) && (zAngle == 0.0);
          vertical = (xAngle == 0.0) && (zAngle == Math.PI/2);
          wide = (yAngle == Math.PI/2) && (zAngle == Math.PI/2);
          if(horizontal) alignment = Constants.HORIZONTAL;
          if(vertical) alignment = Constants.VERTICAL;
          if(wide) alignment = Constants.WIDTH;
        }
        
        public void setXAngle(double t) {
            double pi = Math.PI;
            xAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosX = Math.cos(xAngle);
            sinX = Math.sin(xAngle);
            tanX = Math.tan(xAngle);
            checkAlignment();
        }
        public void setYAngle(double t) {
            double pi = Math.PI;
            yAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosY = Math.cos(yAngle);
            sinY = Math.sin(yAngle);
            tanY = Math.tan(yAngle);
            checkAlignment();
        }
        public void setZAngle(double t) {
            double pi = Math.PI;
            zAngle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            cosZ = Math.cos(zAngle);
            sinZ = Math.sin(zAngle);
            tanZ = Math.tan(zAngle);
            checkAlignment();
        }
        public final double getXAngle() {return(xAngle);}
        public final double getYAngle() {return(yAngle);}
        public final double getZAngle() {return(zAngle);}
        public final double getSinX() {return(sinX);}
        public final double getSinY() {return(sinY);}
        public final double getSinZ() {return(sinZ);}
        public final double getCosX() {return(cosX);}
        public final double getCosY() {return(cosY);}
        public final double getCosZ() {return(cosZ);}
        public final double getTanX() {return(tanX);}
        public final double getTanY() {return(tanY);}
        public final double getTanZ() {return(tanZ);}
        
        public final void setAlignment(Constants.Alignment a) {
            alignment = a;
            if(a == Constants.VERTICAL) {
              setXAngle(0.0);
              setYAngle(0.0);
              setZAngle(Math.PI/2);
            } else if(a == Constants.HORIZONTAL) {
              setXAngle(0.0);
              setYAngle(0.0);
              setZAngle(0.0);
            } else {
              setXAngle(Math.PI/2);
              setYAngle(0.0);
              setZAngle(0.0);
            }
        }
        public final Constants.Alignment getAlignment() {return alignment;}
            
        public final boolean isVertical() {return vertical;}
        public final boolean isHorizontal() {return horizontal;}
        public final boolean isWide() {return wide;}
        
        public final int getThickness() {return thickness;}
        public final void setThickness(int thickness) {this.thickness = thickness;}
        
        public final double getLength() {return length;}
        public final void setLength(double length) {
            this.length = length; 
        }
        
        public final void setLongWall(boolean b) {longWall = b;}
        public final boolean isLongWall() {return longWall;}
        
        public final double getTemperature() {return temperature;}
//        public final void setTemperature(int t) {setTemperature((double)t);}  //for connection to sliders, etc.
        public final void setTemperature(double t) {temperature = t;}
        
        public final boolean isAdiabatic() {return adiabatic;}
        public final void setAdiabatic(boolean a) {adiabatic = a;}
        
        //Meter that evaluates pressure based on accumulated momentum from collisions with wall
        public class MeterPressure extends etomica.Meter {
            
            private double timeSum;
            public MeterPressure() {
                this(Simulation.instance);
            }
            public MeterPressure(Simulation sim) {
                super(sim);
            }
            public etomica.units.Dimension getDimension() {return etomica.units.Dimension.PRESSURE;}

            public void intervalAction(Integrator.IntervalEvent evt) {
                IntegratorMD integrator = (IntegratorMD)evt.getSource();
                if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
                timeSum += integrator.timeStep * integrator.interval;
	            if(--iieCount == 0) {
	                iieCount = updateInterval;
	                updateSums();
	                timeSum = 0.0;
	                pAccumulator = 0.0;
	            }
            }

            public double currentValue() {
                double flux = pAccumulator/timeSum;   //divide by time interval
                flux /= length; //divide by area
                return flux;
            }
        }//end of MeterPressure
    }
    
    //prototype of a real atom type
/*    public final static class Carbon extends Sphere {
        public Carbon() {
            super(12.0, Color.black, 1.1);  //mass, color, diameter  
        }
    }
    public final static class Carbon12 extends Sphere {
        public Carbon12() {
            super(12.0, Color.black, 1.1);
            this.setName("Carbon" + Integer.toString(CarbonID++));
        }
    }
    
    public final static class Hydrogen extends Sphere {
        public Hydrogen() {
            super(1.0, Color.cyan, 0.5);
            this.setName("Hydrogen" + Integer.toString(HydrogenID++));
        }
    }
    
    public final static class Oxygen extends Sphere {
        public Oxygen() {
            super(16.0, Color.red, 1.3);
            this.setName("Oxygen" + Integer.toString(OxygenID++));
        }
    }
    
    public final static class Nitrogen extends Sphere {
        public Nitrogen() {
            super(14.0, Color.blue, 1.2);
            this.setName("Nitrogen" + Integer.toString(NitrogenID++));
        }
    }    
    
*/    
    //interfaces for anisotropic atom types
    public interface Rotator {
        public double[] momentOfInertia(); //diagonal elements of (diagonalized) moment of inertia; should always be a 3-element array
    }
    public interface SphericalTop extends Rotator {} //guarantees Ixx = Iyy = Izz
    public interface CylindricalTop extends Rotator {} //guarantees Ixx = Iyy
    public interface AsymmetricTop extends Rotator {} //all moment-of-inertia elements unequal
    
    /**
     * Type for an AtomGroup atom.
     */
    public static class Group extends AtomType {
        public Group(AtomFactory creator) {
            super(creator);
        }
    }

    private static class Null extends AtomType.Group {
        public Null() {super(null);}
    }
    public static final AtomType.Null NULL = new Null();
        
}
        
