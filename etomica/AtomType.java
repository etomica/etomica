package simulate;
import java.awt.Graphics;
import java.awt.Color;

public interface AtomType {
            
    public void draw(Graphics g, int origin[], double s, Color color, PhaseSpace.AtomCoordinate c);

            
    public final class Disk implements AtomType {
        private double diameter, radius;
        private double mass, rm;
        
        public Disk(double d, double m) {
            setDiameter(d);
            setMass(m);
        }
                    
        public final double getDiameter() {return diameter;}
        public final double getRadius() {return radius;}
        
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
        public final double rm() {return rm;}
        public final void setRm(double rm) {
            this.rm = rm;
            mass = (rm==0.0) ? Double.MAX_VALUE : 1.0/rm;
        }
        
        public final double mass() {return mass;}
        /**
        * Sets diameter of this atom and updates radius accordingly.
        *
        * @param d   new value for diameter
        */
        public final void setDiameter(double d) {diameter = d; radius = 0.5*d;}
        /**
        * Sets radius of this atom and updates diameter accordingly.
        *
        * @param r   new value for radius
        */
        public final void setRadius(double r) {this.setDiameter(2.0*r);}
                
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
        public void draw(Graphics g, int[] origin, double scale, Color color, PhaseSpace.AtomCoordinate c) {
            PhaseSpace2D.AtomCoordinate c2 = (PhaseSpace2D.AtomCoordinate)c;
            double toPixels = scale*DisplayConfiguration.SIM2PIXELS;
            int sigmaP = (int)(toPixels*diameter);
            int xP = origin[0] + (int)(toPixels*(c2.r.x-radius));
            int yP = origin[1] + (int)(toPixels*(c2.r.y-radius));
            g.setColor(color);
            g.fillOval(xP,yP,sigmaP,sigmaP);
        }
    }
}
        
