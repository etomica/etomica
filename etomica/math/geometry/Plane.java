
//(maybe should put this in a Space package hierarchy)
package etomica.math.geometry;

import etomica.Space3D;
import etomica.Space;

/**
 * Class describing a plane occupying a 3-dimensional space.
 *
 * @author David Kofke
 */
 
 //reference:  http://mathworld.wolfram.com/Plane.html
 
 /* History of changes
  * 09/07/02 (DAK) new
  * 09/21/02 (DAK) added method for computing in-plane vectors
  */

public class Plane {
    
    //coefficients defining plane, ax + by + cz + d = 0
    //such that a^2 + b^c + c^2 = 1
    private double a, b, c, d;
    private Space3D.Vector[] inPlane; //work vectors used by inPlaneSquare method
    private Space3D.Vector ctr; //work vector used by inPlaneSquare method (no-x0 version)
    
    /**
     * Tolerance used to judge if a given point is in the plane.
     * Method inPlane is true if the point given to is at a distance
     * from the plane less than this value.  Default is 1e-5.
     */
    public double epsilon = 1.0e-5;
    
    /**
     * Default constructor returns the y-z plane.
     */
    public Plane() {
        this(1.0, 0.0, 0.0, 0.0);
    }
    /**
     * Constructs a new Plane as a copy of the given one.
     */
    public Plane(Plane original) {
        this(original.a, original.b, original.c, original.d);
    }
    
    /**
     * Constructs a Plane satisfying the equation a x + b y + c z + d = 0
     */
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
     * Sets the given vector to be a unit normal vector to the plane and returns it.  
     * If argument is null, creates a new 3D Vector to return.
     */
    public Space3D.Vector getNormalVector(Space3D.Vector normal) {
        if(normal == null) return new Space3D.Vector(a, b, c);
        else {
            normal.E(a,b,c);
            return normal;
        }
    }
    /**
     * Computes and returns two unit vectors perpendicular to each other and in the plane.
     * Sets vector in given array to the vectors if it is not null and is of length 2;
     * otherwise creates new array of in-plane vectors and returns it.  If plane is
     * parallel to a coordinate-axes plane, returns coordinate vectors (e.g., if parallel
     * to x-y plane, returns unit vectors in x and y, respectively).  Otherwise, returns
     * one vector that follows the intersection of this plane and the x-y plane; the
     * second vector is (uniquely) defined as perpendicular to both the normal and the first
     * in-plane vector.
     */
    public Space3D.Vector[] inPlaneVectors(Space3D.Vector[] p) {
        if(p == null || p.length != 2) p = new Space3D.Vector[] {new Space3D.Vector(), new Space3D.Vector()};
        Space3D.Vector p1 = p[0];
        Space3D.Vector p2 = p[1];
        if(a == 0 && b == 0) {//parallel to xy plane
            p1.E(1, 0, 0);
            p2.E(0, 1, 0);
        } else {
            p1.E(b, -a, 0);
            p1.normalize();
            if(c == 0) p2.E(0, 0, 1);
            else {
                p2.E(a, b, -(a*a+b*b)/c);
                p2.normalize();
            }
        }
        return p;
    }
    
    /**
     * Sets the orientation of the plane to be normal to the given vector.
     */
    public void setNormalVector(Space.Vector n) {
        if(!(n instanceof Space3D.Vector)) throw new IllegalArgumentException("Error: Plane.setNormalVector requires 3D vector as argument");
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
     * Sets the given vector to be the point in the plane closest to the origin.
     * If vector is null, makes a new Vector and returns it.
     */
    public Space3D.Vector center(Space3D.Vector v) {
        if(v == null) return new Space3D.Vector(-d*a, -d*b, -d*c);
        else {
            v.E(-d*a, -d*b, -d*c);
            return v;
        }
    }
    
    /**
     * Determines the point in the plane closest to the point x0 (first argument).
     * Uses second vector to hold result that is returned; if second vector is 
     * null, creates a new Vector instance.
     */
    public Space3D.Vector nearestPoint(Space3D.Vector x0, Space3D.Vector point) {
        if(point == null) point = new Space3D.Vector();
        double factor = distanceTo(x0);
        point.E(x0.x(0)-factor*a, x0.x(1)-factor*b, x0.x(2)-factor*c);
        return point;
    }

    /**
     * Perpendicular distance from the plane to the given point.
     */
    public double distanceTo(Space3D.Vector x0) {
        return a*x0.x(0) + b*x0.x(1) + c*x0.x(2) + d;
    }
    
    /**
     * Shifts the plane at fixed orientation so that it contains the given point.
     */
    public void moveTo(Space3D.Vector r) {
        d = -(a*r.x(0) + b*r.x(1) + c*r.x(2));
    }        
    
    /**
     * Dihedral angle between this plane and the given plane.
     */
    public double dihedralAngle(Plane p) {
        Space3D.Vector n = p.getNormalVector(null);
        return a*n.x(0) + b*n.x(1) + c*n.x(2);
    }
    
    /**
     * Sets the parameters of this plane to make it equivalent to the given one.
     */
    public void E(Plane p) {
        a = p.a;
        b = p.b;
        c = p.c;
        d = p.d;
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
     * Returns true if the given point is within a distance epsilon of this plane.
     */
    public boolean inPlane(Space3D.Vector p) {
        return Math.abs(distanceTo(p)) < epsilon;
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
    
    /**
     * Returns four points as the vertices of a square of size d centered on x0.  Uses 
     * and returns the given array of vectors if it is not null and is of length 4; 
     * otherwise makes a new array and returns it.  Square is aligned so that its vertices
     * fall on the lines defined by the inPlaneVectors result.
     */
    public Space3D.Vector[] inPlaneSquare(Space3D.Vector x0, double d, Space3D.Vector[] s) {
        if(s == null || s.length != 4) s = new Space3D.Vector[] {new Space3D.Vector(), new Space3D.Vector(), new Space3D.Vector(), new Space3D.Vector()};
        inPlane = inPlaneVectors(inPlane);
        Space3D.Vector p1 = inPlane[0];
        Space3D.Vector p2 = inPlane[1];
      //  System.out.println("ctr: "+x0.toString());
      //  System.out.println("p1 : "+p1.toString());
      //  System.out.println("p2 : "+p2.toString());
        double t = d/Math.sqrt(2.0);
        s[0].E(x0); s[0].PEa1Tv1(+t, p1);
        s[1].E(x0); s[1].PEa1Tv1(-t, p1);
        s[2].E(x0); s[2].PEa1Tv1(+t, p2);
        s[3].E(x0); s[3].PEa1Tv1(-t, p2);
        return s;
    }
    
    /**
     * inPlaneSquare with square centered on point in plane closest to origin.
     */
    public Space3D.Vector[] inPlaneSquare(double d, Space3D.Vector[] s) {
        return inPlaneSquare(center(ctr), d, s);
    }
    
     
    //resets coefficients so that a^2 + b^2 + c^2 = 1;
    private void normalize() {
        double norm = Math.sqrt(a*a + b*b + c*c);
        a /= norm;
        b /= norm;
        c /= norm;
        d /= norm;
    }
    
    public static void main(String[] args) {
        Plane plane = new Plane(2.8, 7.5, -2.1, 293.8);
    //    Plane plane = new Plane(5.0, 0.0, 0.0, 293.8);
        Space3D.Vector[] p = plane.inPlaneVectors(null);
        Space3D.Vector n = plane.getNormalVector(null);
        System.out.println("Normal: "+n.toString());
        System.out.println("plane1: "+p[0].toString());
        System.out.println("plane2: "+p[1].toString());
        System.out.println("p1.n: "+p[0].dot(n));
        System.out.println("p2.n: "+p[1].dot(n));
        System.out.println("p1.p2: "+p[0].dot(p[1]));
    }
    
}//end of Plane