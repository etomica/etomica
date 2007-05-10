
//(maybe should put this in a Space package hierarchy)
package etomica.math.geometry;

import etomica.space.IVector;
import etomica.space3d.IVector3D;
import etomica.space3d.Vector3D;

/**
 * Class describing a plane occupying a 3-dimensional space.  
 *
 * @author David Kofke
 */
 
 //reference:  http://mathworld.wolfram.com/Plane.html

public class Plane implements java.io.Serializable {
    
    private static final long serialVersionUID = 1L;
    //coefficients defining plane, ax + by + cz + d = 0
    //such that a^2 + b^c + c^2 = 1
    private double a, b, c, d;
    private IVector3D[] inPlane; //work vectors used by inPlaneSquare method
    private IVector3D work0; //work vector used by inPlaneSquare method (no-x0 version)
    private IVector3D work1; //work vector used by inPlaneSquare method (no-x0 version)
    
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
     * Defines the plane via a normal vector (1st argument) and
     * a point (2nd argument).
     */
    public void setNormalPoint(IVector3D normal, IVector3D point) {
        a = normal.x(0);
        b = normal.x(1);
        c = normal.x(2);
        d = -normal.dot(point);
        normalize();
    }
    
    /**
     * Defines the plane by specifying three points that lie in it.
     */
    public void setThreePoints(Vector3D p1, Vector3D p2, Vector3D p3) {
        if(work0 == null) {
            work0 = new Vector3D();
            work1 = new Vector3D();
        }
        work0.Ev1Mv2(p2,p1);
        work1.Ev1Mv2(p3,p1);
        work1.XE(work0);
        work1.normalize();
        setNormalPoint(work1, p1);
//        work.E(1.0);
//        Vector3D px = new Vector3D(p1.x(0),p2.x(0),p3.x(0));
//        Vector3D py = new Vector3D(p1.x(1),p2.x(1),p3.x(1));
//        Vector3D pz = new Vector3D(p1.x(2),p2.x(2),p3.x(2));
//        double aa = det(work, py, pz);
//        double bb = det(px, work, pz);
//        double cc = det(px, py, work);
//        double dd = -det(px, py, pz);
//        double norm = Math.sqrt(aa*aa+bb*bb+cc*cc);
//        aa/=norm; bb/=norm; cc/=norm; dd/=norm;
//        a = aa; b = bb; c = cc; d = dd;
//        System.out.println("a "+a+" "+aa);
//        System.out.println("b "+b+" "+bb);
//        System.out.println("c "+c+" "+cc);
//        System.out.println("d "+d+" "+dd);
//        System.out.println("check1 "+(a*p1.x(0)+b*p1.x(1)+c*p1.x(2)+d));
//        System.out.println("check2 "+(a*p2.x(0)+b*p2.x(1)+c*p2.x(2)+d));
//        System.out.println("check3 "+(a*p3.x(0)+b*p3.x(1)+c*p3.x(2)+d));
//        normalize();
    }
    
    //used by setThreePoints
//    private double det(Vector3D p1, Vector3D p2, Vector3D p3) {
//        return p1.x(0)*(p2.x(1)*p3.x(2)-p3.x(1)*p2.x(2))
//             - p1.x(1)*(p2.x(0)*p3.x(2)-p3.x(0)*p2.x(2))
//             + p1.x(2)*(p2.x(0)*p3.x(1)-p3.x(0)*p2.x(1));
//    }
    
    
    /**
     * Sets the given vector to be a unit normal vector to the plane and returns it.  
     * If argument is null, creates a new 3D Vector to return.
     */
    public IVector getNormalVector(IVector3D normal) {
        if(normal == null) return new Vector3D(a, b, c);
        normal.E(a,b,c);
        return normal;
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
    public IVector3D[] inPlaneVectors(IVector3D[] p) {
        if(p == null || p.length != 2) p = new Vector3D[] {new Vector3D(), new Vector3D()};
        IVector3D p1 = p[0];
        IVector3D p2 = p[1];
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
    public void setNormalVector(IVector3D n) {
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
    public IVector3D center(IVector3D v) {
        if(v == null) return new Vector3D(-d*a, -d*b, -d*c);
        v.E(-d*a, -d*b, -d*c);
        return v;
    }
    
    /**
     * Determines the point in the plane closest to the point x0 (first argument).
     * Uses second vector to hold result that is returned; if second vector is 
     * null, creates a new Vector instance.
     */
    public IVector3D nearestPoint(Vector3D x0, Vector3D point) {
        if(point == null) point = new Vector3D();
        double factor = distanceTo(x0);
        point.E(x0.x(0)-factor*a, x0.x(1)-factor*b, x0.x(2)-factor*c);
        return point;
    }

    /**
     * Perpendicular distance from the plane to the given point.
     */
    public double distanceTo(Vector3D x0) {
        return a*x0.x(0) + b*x0.x(1) + c*x0.x(2) + d;
    }
    
    /**
     * Shifts the plane at fixed orientation so that it contains the given point.
     */
    public void moveTo(Vector3D r) {
        d = -(a*r.x(0) + b*r.x(1) + c*r.x(2));
    }        
    
    /**
     * Dihedral angle between this plane and the given plane.
     */
    public double dihedralAngle(Plane p) {
        IVector n = p.getNormalVector(null);
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
     * plane or over the side of the plane by epsilon or less
     * toward which the normal vector points.  The direction
     * of the normal vector can be inverted using the invert
     * method.
     */
    public boolean isPositiveSide(Vector3D p) {
        return distanceTo(p) > epsilon;
    }
    
    /**
     * Returns true if the given point is within a distance epsilon of this plane.
     */
    public boolean inPlane(Vector3D p) {
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
    public IVector3D[] inPlaneSquare(IVector x0, double d, IVector3D[] s) {
        if(s == null || s.length != 4) s = new IVector3D[] {new Vector3D(), new Vector3D(), new Vector3D(), new Vector3D()};
        inPlane = inPlaneVectors(inPlane);
        IVector p1 = inPlane[0];
        IVector p2 = inPlane[1];
      //  System.out.println("work: "+x0.toString());
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
    public IVector3D[] inPlaneSquare(double d, IVector3D[] s) {
        return inPlaneSquare(center(work0), d, s);
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
        IVector[] p = plane.inPlaneVectors(null);
        IVector n = plane.getNormalVector(null);
        System.out.println("Normal: "+n.toString());
        System.out.println("plane1: "+p[0].toString());
        System.out.println("plane2: "+p[1].toString());
        System.out.println("p1.n: "+p[0].dot(n));
        System.out.println("p2.n: "+p[1].dot(n));
        System.out.println("p1.p2: "+p[0].dot(p[1]));
    }
    
}//end of Plane
