package etomica;
import java.awt.Graphics;
import java.awt.Color;
import etomica.units.*;
import etomica.electrostatics.*;

/**
 * AtomType holds atom parameters, including the mass.  
 * It is used to set the general features of the atom (e.g., whether it is a disk, wall, sphere, etc.),
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
    public static String getVersion() {return "01.03.05.0";}
    private double mass, rm;
    private Color color;
    private ElectroType electroType;
    private static int CarbonID = 0;
    private static int HydrogenID = 0;
    private static int OxygenID = 0;
    private static int NitrogenID = 0;
    private String name;
    
    public AtomType(double m, Color c) {
        setMass(m);
        setColor(c);
        setElectroType(new ElectroType.Null());
    }
    
    /**
     * Returns default coordinate type, which has no orientational component.
     * Override for atom types that require other coordinate features.
     */
    public Space.Coordinate makeCoordinate(Atom a) {
        return a.parentSimulation().space().makeCoordinate(a);
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

    // AtomType color may or may not determine atom.color (and thus drawn color), depending on ColorScheme
    public final Color color() {return color;}
    public final void setColor(Color c) {color = c;}
    
    public abstract void draw(Graphics g, int origin[], double t, Atom atom);

    /**
     * Accessor method of the name of this species
     * 
     * @return The given name of this species
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this species
     * The species' name provides a convenient way to label output data that is 
     * associated with this species.  This method might be used, for example, to place
     * a heading on a column of data.
     * Default name is "Species" followed by the integer species index of this species.
     * 
     * @param name The name string to be associated with this species
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the species
     */
    public String toString() {return getName();}  //override Object method

    // Rod-shaped atom.  Assumed to be in 1D
    public static class Rod extends Disk {
        
        private int dh2;
        public Rod(double m, Color c, double d) {
            super(m, c, d);
            dh2 = Space1D.drawingHeight/2;
        }
        
        /**
        * Draws this atom using current values of its position, diameter and color.
        * Drawing position is determined as follows.  The atoms coordinates in
        * Angstroms are converted to pixels by applying a scaling factor; these
        * drawing coordinates may be shifted by some amount as given by the array
        * <code>origin</code> before the atom is drawn.
        *
        * @param g         graphics object to which atom is drawn
        * @param origin    origin of drawing coordinates (pixels)
        * @param scale     factor determining size of drawn image relative to
        *                  nominal drawing size
        */
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {
            Space.Vector r = atom.coordinate().position();
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(r.component(0)-radius));
            int yP = origin[1]-dh2;
            g.setColor(atom.color);
            g.fillRect(xP,yP,sigmaP,Space1D.drawingHeight);
            electroType.draw(g, origin, toPixels, r);
        }
    }

    // Disk-shaped atom.  Assumed to be in 2D
    public static class Disk extends AtomType {
        
        double diameter, radius;
        
        public Disk(double m, Color c, double d) {
            super(m, c);
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
                
        /**
        * Draws this atom using current values of its position, diameter and color.
        * Drawing position is determined as follows.  The atoms coordinates in
        * Angstroms are converted to pixels by applying a scaling factor; these
        * drawing coordinates may be shifted by some amount as given by the array
        * <code>origin</code> before the atom is drawn.
        *
        * @param g         graphics object to which atom is drawn
        * @param origin    origin of drawing coordinates (pixels)
        * @param scale     factor determining size of drawn image relative to
        *                  nominal drawing size
        */
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {
            Space.Vector r = atom.coordinate().position();
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(r.component(0)-radius));
            int yP = origin[1] + (int)(toPixels*(r.component(1)-radius));
            g.setColor(atom.color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
            electroType.draw(g, origin, toPixels, r);
        }
    }
    
    /**
     * Atom type for a sphere that has some feature depending upon an orientation coordinate.
     * For example an orientational dependent potential may be attached to an otherwise spherical atom
     */
    public final static class OrientedSphere extends Disk implements SphericalTop {
        
        private final double[] I = new double[3];
        public OrientedSphere(double m, Color c, double d) {
            super(m,c,d);
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
        
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {
            Space.Vector r = atom.coordinate().position();
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(r.component(0)-radius));
            int yP = origin[1] + (int)(toPixels*(r.component(1)-radius));
            g.setColor(atom.color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
            g.setColor(Color.red);
            double theta = ((Space.Coordinate.Angular)atom.coordinate()).orientation().angle()[0];
            int dx = (int)(toPixels*radius*Math.cos(theta));
            int dy = (int)(toPixels*radius*Math.sin(theta));
            int dxy = (int)(toPixels*radius);
            xP += dxy; yP += dxy;
            g.drawLine(xP-dx, yP-dy, xP+dx, yP+dy);
            
            electroType.draw(g, origin, toPixels, r);
        }
    }
    
    // Disk with a concentric well.  Assumed to be in 2D
    public final static class Well extends Disk {  
        
        private double lambda;                    //diameter of well, in units of core diameter
        private double wellDiameter, wellRadius;  //size of well, in simulation units
        
        public Well(double m, Color c, double d, double l) {
            super(m, c, d);
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
                
        /**
        * Draws this atom using current values of its position, diameter and color.
        * Drawing position is determined as follows.  The atoms coordinates in
        * Angstroms are converted to pixels by applying a scaling factor; these
        * drawing coordinates may be shifted by some amount as given by the array
        * <code>origin</code> before the atom is drawn.
        *
        * @param g         graphics object to which atom is drawn
        * @param origin    origin of drawing coordinates (pixels)
        * @param scale     factor determining size of drawn image relative to
        *                  nominal drawing size
        */
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {
            Space.Vector r = atom.coordinate().position();
 //           int yP = 0;
            //Draw core
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(r.component(0)-radius));
            int yP = origin[1] + (int)(toPixels*(r.component(1)-radius));  //removed for 1D
            g.setColor(atom.color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
            
            //Draw well
            sigmaP = (int)(toPixels*wellDiameter);
            xP = origin[0] + (int)(toPixels*(r.component(0)-wellRadius));
            yP = origin[1] + (int)(toPixels*(r.component(1)-wellRadius));  //removed for 1D
            g.setColor(Color.lightGray);
            g.drawOval(xP,yP,sigmaP,sigmaP);
            electroType.draw(g, origin, toPixels, r);
        }
    }

    // Wall-shaped atom.  Assumed to be in Space2D
    public final static class Wall extends AtomType {
        
        int thickness = 4;  //thickness when drawn to screen (if horizontal or vertical)
        private boolean vertical, horizontal;
        private double cosA, sinA;
//        double[] f = new double[Space.D];   //force on wall
//        double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
        protected double angle;   //orientation (radians) with respect to positive x-axis, taking r at the origin
        protected double length;  //length of line in simulation units
        protected boolean longWall;  //set to true if wall is infinite in length
        //  these may belong in integratorAgent
        protected double temperature = Kelvin.UNIT.toSim(300.);
        protected boolean adiabatic = true;
        private Constants.Alignment alignment;
        
        public Wall(double m, Color c, double l, double a) {
            super(m, c);
            setLength(l);
            setAngle(a);
        }
                    
        public void setAngle(double t) {
            double pi = Math.PI;
            angle = (t <= 2.*pi) ? t : (t % (2.0*pi));
            horizontal = (angle == 0.0) || (Math.abs(angle) == pi);
            vertical = (Math.abs(angle) == pi/2.) || (Math.abs(angle) == 1.5*pi);
            if(horizontal) alignment = Constants.HORIZONTAL;
            if(vertical) alignment = Constants.VERTICAL;
            cosA = Math.cos(angle);
            sinA = Math.sin(angle);
        }
        public final double getAngle() {return angle;}
        
        public final void setAlignment(Constants.Alignment a) {
            alignment = a;
            if(a == Constants.VERTICAL) {setAngle(Math.PI/2.);}
            else if(a == Constants.HORIZONTAL) {setAngle(0.0);}
        }
        public final Constants.Alignment getAlignment() {return alignment;}
            
        public final boolean isVertical() {return vertical;}
        public final boolean isHorizontal() {return horizontal;}
        
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
     
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {
            Space.Vector r = atom.coordinate().position();
            int xP = origin[0] + (int)(toPixels*r.component(0));
            int yP = origin[1] + (int)(toPixels*r.component(1));
            int t = Math.max(1,(int)((double)thickness*(double)toPixels/(double)BaseUnit.Length.Sim.TO_PIXELS));
            g.setColor(atom.color);
            if(!(horizontal || vertical)) {  //not horizontal or vertical; draw line
                int x1 = xP + (int)(toPixels*length*cosA);
                int y1 = yP + (int)(toPixels*length*sinA);
                g.drawLine(xP, yP, x1, y1);
            }
            else if(((Wall)atom.type).isLongWall()) {
                java.awt.Rectangle rect = g.getClipBounds();
          //      int wP = vertical ? t : (int)(toPixels*atom.parentPhase().boundary().dimensions().component(1));
          //      int hP = horizontal ? t : (int)(toPixels*atom.parentPhase().boundary().dimensions().component(0));
                int wP = vertical ? t : Integer.MAX_VALUE;
                int hP = horizontal ? t : Integer.MAX_VALUE;
           //     int X = vertical ? xP : origin[0];
           //     int Y = horizontal ? yP : origin[1];
                int X = vertical ? xP : 0;
                int Y = horizontal ? yP : 0;
                g.fillRect(X,Y,wP,hP);
            }   
            else {                           //horizontal or vertical; draw box
                int wP = vertical ? t : (int)(toPixels*length);
                int hP = horizontal ? t : (int)(toPixels*length);
                g.fillRect(xP,yP,wP,hP);
            }
        }
    }
    
    //prototype of 3D type
    public static class Sphere extends AtomType {
        
        private double diameter, radius;
        
        public Sphere(double m, Color c, double d) {
            super(m, c);
            setDiameter(d);
        }
        
        public final double diameter() {return diameter;}
        public final double radius() {return radius;}
        public final void setDiameter(double d) {diameter = d; radius = 0.5*d;}
        public final void setRadius(double r) {this.setDiameter(2.0*r);}
        
        //Need to complete this
        public void draw(Graphics g, int[] origin, double toPixels, Atom atom) {}
    }
    
    //prototype of a real atom type
    public final static class Carbon extends Sphere {
        public Carbon() {
            super(12.0, Color.black, 1.1);  //mass, color, diameter  
        }
    }
    public final static class Carbon12 extends Disk {
        public Carbon12() {
            super(12.0, Color.black, 1.1);
            this.setName("Carbon" + Integer.toString(CarbonID++));
        }
    }
    
    public final static class Hydrogen extends Disk {
        public Hydrogen() {
            super(1.0, Color.cyan, 0.5);
            this.setName("Hydrogen" + Integer.toString(HydrogenID++));
        }
    }
    
    public final static class Oxygen extends Disk {
        public Oxygen() {
            super(16.0, Color.red, 1.3);
            this.setName("Oxygen" + Integer.toString(OxygenID++));
        }
    }
    
    public final static class Nitrogen extends Disk {
        public Nitrogen() {
            super(14.0, Color.blue, 1.2);
            this.setName("Nitrogen" + Integer.toString(NitrogenID++));
        }
    }    
    
    
    //interfaces for anisotropic atom types
    public interface Rotator {
        public double[] momentOfInertia(); //diagonal elements of (diagonalized) moment of inertia; should always be a 3-element array
    }
    public interface SphericalTop extends Rotator {} //guarantees Ixx = Iyy = Izz
    public interface CylindricalTop extends Rotator {} //guarantees Ixx = Iyy
    public interface AsymmetricTop extends Rotator {} //all moment-of-inertia elements unequal
        
}
        
