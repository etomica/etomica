package simulate;
import java.awt.Graphics;
import java.awt.Color;

/**
 * AtomType holds atom parameters, including the mass.  
 * It is used to set the general features of the atom (e.g., whether it is a disk, wall, sphere, etc.),
 * and to prescribe how it is drawn to the screen.
 * The AtomType of an atom is set by Species when it constructs a molecule.  Each Atom has an instance variable
 * named "type that holds the AtomType object; this may be different for each atom in a molecule, or
 * it may refer to a common AtomType object, as prescribed by the Species.
 * AtomType could also be used to define particular elemental atoms (Carbon, Oxygen, etc.).
 * 
 */

public abstract class AtomType {
    private double mass, rm;
    private Color color;
    
    public AtomType(double m, Color c) {
        setMass(m);
        setColor(c);
    }
            
    /**
    * Sets  mass of this atom and updates reciprocal mass accordingly.  Setting
    * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass 
    * to be set to zero.
    * 
    * @param mass   new value for mass
    */
    public final void setMass(double mass) {
        this.mass = mass; 
        rm = (mass==Double.MAX_VALUE) ? 0.0 : 1.0/mass;
    }
    public final void setRm(double rm) {
        this.rm = rm;
        mass = (rm==0.0) ? Double.MAX_VALUE : 1.0/rm;
    }
    public final double rm() {return rm;}
    public final double mass() {return mass;}
    
    // AtomType color may or may not determine atom.color (and thus drawn color), depending on ColorScheme in Species
    public final Color color() {return color;}
    public final void setColor(Color c) {color = c;}
    
    public abstract void draw(Graphics g, int origin[], double s, Color color, Space.AtomCoordinate c);

    
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
        public void draw(Graphics g, int[] origin, double scale, Color color, Space.AtomCoordinate c) {
            Space2D.AtomCoordinate c2 = (Space2D.AtomCoordinate)c;
            double toPixels = scale*DisplayConfiguration.SIM2PIXELS;
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(c2.r.x-radius));
            int yP = origin[1] + (int)(toPixels*(c2.r.y-radius));
            g.setColor(color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
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
        
        public final void setDiameter(double d) {super.setDiameter(d); setLambda(lambda);}
        public final void setLambda(double l) {lambda = l; wellDiameter = lambda*diameter(); wellRadius = 0.5*wellDiameter;}
                
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
        public void draw(Graphics g, int[] origin, double scale, Color color, Space.AtomCoordinate c) {
            Space2D.AtomCoordinate c2 = (Space2D.AtomCoordinate)c;
            double toPixels = scale*DisplayConfiguration.SIM2PIXELS;

            //Draw core
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(c2.r.x-radius));
            int yP = origin[1] + (int)(toPixels*(c2.r.y-radius));
            g.setColor(color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
            
            //Draw well
            sigmaP = (int)(toPixels*wellDiameter);
            xP = origin[0] + (int)(toPixels*(c2.r.x-wellRadius));
            yP = origin[1] + (int)(toPixels*(c2.r.y-wellRadius));
            g.setColor(color);
            g.drawOval(xP,yP,sigmaP,sigmaP);
        }
    }

    // Wall-shaped atom.  Assumed to be in Space2D
    public final static class Wall extends AtomType {
        
        int thickness = 4;  //thickness when drawn to screen (if horizontal or vertical)
        private boolean vertical, horizontal;
        private double cosA, sinA;
//        double[] f = new double[Space.D];   //force on wall
//        double[] r1 = new double[Space.D];  //other end of line in simulation units (first end is denoted r)
        protected int angle;   //orientation (degrees) with respect to positive x-axis, taking r at the origin
        protected double length;  //length of line in simulation units
        protected boolean longWall;  //set to true if wall is infinite in length
        //  these may belong in integratorAgent
        protected double temperature = 300.;
        protected boolean adiabatic = true;
        
        public Wall(double m, Color c, double l, int a) {
            super(m, c);
            setLength(l);
            setAngle(a);
        }
                    
        public void setAngle(int t) {
            angle = (t <= 360) ? t : (t % 360);
            horizontal = (angle == 0) || (Math.abs(angle) == 180);
            vertical = (Math.abs(angle) == 90) || (Math.abs(angle) == 270);
            cosA = Math.cos((double)angle*Math.PI/180.);
            sinA = Math.sin((double)angle*Math.PI/180.);
        }
        public final int getAngle() {return angle;}
            
        public final boolean isVertical() {return vertical;}
        public final boolean isHorizontal() {return horizontal;}
        
        public final int getThickness() {return thickness;}
        public final void setThickness(int thickness) {this.thickness = thickness;}
        
        public final double getLength() {return length;}
        public final void setLength(double length) {this.length = length; longWall = (length == Double.MAX_VALUE);}
        
        public final boolean isLongWall() {return longWall;}
        
        public final double getTemperature() {return temperature;}
        public final void setTemperature(int t) {setTemperature((double)t);}  //for connection to sliders, etc.
        public final void setTemperature(double t) {temperature = t;}
        
        public final boolean isAdiabatic() {return adiabatic;}
        public final void setAdiabatic(boolean a) {adiabatic = a;}
     
        public void draw(Graphics g, int[] origin, double scale, Color color, Space.AtomCoordinate c) {
            Space2D.AtomCoordinate c2D = (Space2D.AtomCoordinate)c;
            double toPixels = scale*DisplayConfiguration.SIM2PIXELS;
            int xP = origin[0] + (int)(toPixels*c2D.r.x);
            int yP = origin[1] + (int)(toPixels*c2D.r.y);
            g.setColor(color);
            if(!(horizontal || vertical)) {  //not horizontal or vertical; draw line
                int x1 = xP + (int)(toPixels*length*cosA);
                int y1 = yP + (int)(toPixels*length*sinA);
                g.drawLine(xP, yP, x1, y1);
            }
            else {                           //horizontal or vertical; draw box
                int wP = vertical ? thickness : (int)(toPixels*length);
                int hP = horizontal ? thickness : (int)(toPixels*length);
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
        public void draw(Graphics g, int[] origin, double scale, Color color, Space.AtomCoordinate c) {}
    }
    
    //prototype of a real atom type
    public final static class Carbon extends Sphere {
        public Carbon() {
            super(12.0, Color.black, 1.1);  //mass, color, diameter  
        }
    }
}
        
