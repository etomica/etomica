package etomica.math.geometry;

import etomica.Space3D;

/**
 * Class describing a plane occupying a 3-dimensional space.
 *
 * @author David Kofke
 */
 
 //reference:  http://mathworld.wolfram.com/Plane.html
 
 /* History of changes
  * 09/07/02 (DAK) new
  */

public class Plane {
    
    //coefficients defining plane, ax + by + cz + d = 0
    //such that a^2 + b^c + c^2 = 1
    private double a, b, c, d;
    
    public Plane(double a, double b, double c, double d) {
        if(a == 0 && b == 0 && c == 0) throw new IllegalArgumentException("Arguments to Plane constructor do not define a plane"); 
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
        normalize();
    }
    
    /**
     * Returns a new plane having x, y, and z intercepts given by the arguments.
     */
    public static Plane newInterceptForm(double a, double b, double c) {
        if(a == 0 || b == 0 || c == 0) throw new IllegalArgumentException("Arguments to intercept-form Plane constructor do not define a plane"); 
        return new Plane(1.0/a, 1.0/b, 1.0/c, -1.0);
    }
    
    /**
     * Returns a new plane with the given normal (1st argument)through the 
     * given point (2nd argument).
     */
    public static Plane newNormalPointForm(Space3D.Vector normal, Space3D.Vector point) {
        double a = normal.x(0);
        double b = normal.x(1);
        double c = normal.x(2);
        double x0 = point.x(0);
        double y0 = point.x(1);
        double z0 = point.x(2);
        return new Plane(a, b, c, -(a*x0+b*y0+c*z0));
    }
    
    /**
     * Returns a unit normal vector to the plane.
     */
    public Space3D.Vector getNormalVector() {
        return new Space3D.Vector(a, b, c);
    }
    /**
     * Sets the orientation of the plane to be normal to the given vector.
     */
    public void setNormalVector(Space3D.Vector n) {
        if(n.squared() == 0.0) throw new IllegalArgumentException("Error: attempt to set orientation of plane with respect to an ill-defined vector");
        a = n.x(0);
        b = n.x(1);
        c = n.x(2);
        normalize();
    }
    
    /**
     * Perpendicular distance from the plane to the origin.
     */
    public double getDistanceToOrigin() {return d;}
    
    /**
     * Sets perpendicular distance from the plane to the origin.
     */
    public void setDistanceToOrigin(double d) {this.d = d;}

    /**
     * Reorients and positions the plane so that its intercepts with the
     * x-, y-, and z-axes are the given values, respectively.
     */
    public void setIntercepts(double u, double v, double w) {
        if(u == 0 || v == 0 || w == 0) throw new IllegalArgumentException("Arguments to setIntercept do not define a plane"); 
        a = 1.0/u;
        b = 1.0/v;
        c = 1.0/w;
        d = -1.0;
        normalize();
    }

    /**
     * Returns the intercept of the plane with the x-axis.
     */
    public double xIntercept() {return -d/a;}
    /**
     * Returns the intercept of the plane with the y-axis.
     */
    public double yIntercept() {return -d/b;}
    /**
     * Returns the intercept of the plane with the z-axis.
     */
    public double zIntercept() {return -d/c;}
    
    /**
     * Perpendicular distance from the plane to the given point.
     */
    public double distanceTo(Space3D.Vector x0) {
        return a*x0.x(0) + b*x0.x(1) + c*x0.x(2) + d;
    }
    
    /**
     * Dihedral angle between this plane and the given plane.
     */
    public double dihedralAngle(Plane p) {
        Space3D.Vector n = p.getNormalVector();
        return a*n.x(0) + b*n.x(1) + c*n.x(2);
    }
    
    /**
     * Returns true if the given point is on the side of the 
     * plane toward which the normal vector points.  The direction
     * of the normal vector can be inverted using the invert method.
     */
    public boolean isPositiveSide(Space3D.Vector p) {
        return distanceTo(p) > 0.0;
    }
    
    /**
     * Changes the direction of the normal vector so that it points
     * toward the other side of the plane from its present orientation.
     * Does not affect the location or absolute orientation of the plane.
     */
    public void invert() {
        a = -a;
        b = -b;
        c = -c;
        d = -d;
    }
    
     
    //resets coefficients so that a^2 + b^2 + c^2 = 1;
    private void normalize() {
        double norm = Math.sqrt(a*a + b*b + c*c);
        a /= norm;
        b /= norm;
        c /= norm;
        d /= norm;
    }
    
}//end of Plane