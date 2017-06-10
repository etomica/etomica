/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


//(maybe should put this in a Space package hierarchy)
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Space3D;

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
    private Vector[] inPlane; //work vectors used by inPlaneSquare method
    private Vector work0; //work vector used by inPlaneSquare method (no-x0 version)
    private Vector work1; //work vector used by inPlaneSquare method (no-x0 version)
    private final Space space;
    
    /**
     * Tolerance used to judge if a given point is in the plane.
     * Method inPlane is true if the point given to is at a distance
     * from the plane less than this value.  Default is 1e-5.
     */
    public double epsilon = 1.0e-5;
    
    /**
     * Default constructor returns the y-z plane.
     */
    public Plane(Space space) {
        this(space, 1.0, 0.0, 0.0, 0.0);
    }
    /**
     * Constructs a new Plane as a copy of the given one.
     */
    public Plane(Plane original) {
        this(original.space, original.a, original.b, original.c, original.d);
    }
    
    /**
     * Constructs a Plane satisfying the equation a x + b y + c z + d = 0
     */
    public Plane(Space space, double a, double b, double c, double d) {
        if(a == 0 && b == 0 && c == 0) throw new IllegalArgumentException("Arguments to Plane constructor do not define a plane");
        this.space = space;
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
        normalize();
    }
    
    public double getA() {
        return a;
    }
    
    public double getB() {
        return b;
    }
    
    public double getC() {
        return c;
    }
    
    public double getD() {
        return d;
    }
    
    /**
     * Defines the plane via a normal vector (1st argument) and
     * a point (2nd argument).
     */
    public void setNormalPoint(Vector normal, Vector point) {
        a = normal.getX(0);
        b = normal.getX(1);
        c = normal.getX(2);
        d = -normal.dot(point);
        normalize();
    }
    
    /**
     * Defines the plane by specifying three points that lie in it.
     */
    public void setThreePoints(Vector p1, Vector p2, Vector p3) {
        if(work0 == null) {
            work0 = space.makeVector();
            work1 = space.makeVector();
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
    public void setToNormalVector(Vector normal) {
        normal.setX(0, a);
        normal.setX(1, b);
        normal.setX(2, c);
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
    public void setToInPlaneVectors(Vector[] p) {
        Vector p1 = p[0];
        Vector p2 = p[1];
        work1.E(0);
        if(a == 0 && b == 0) {//parallel to xy plane
            work1.setX(0,1);
            p1.E(work1);
            work1.setX(0,0);
            work1.setX(1,1);
            p2.E(work1);
        } else {
            p1.setX(0, b);
            p1.setX(1, -a);
            p1.setX(2, 0);
            p1.normalize();
            if(c == 0) {
                work1.setX(2,1);
                p2.E(work1);
            }
            else {
                p2.setX(0,a);
                p2.setX(1,b);
                p2.setX(2, -(a*a+b*b)/c);
                p2.normalize();
            }
        }
    }
    
    /**
     * Sets the orientation of the plane to be normal to the given vector.
     */
    public void setNormalVector(Vector n) {
        if(n.isZero()) throw new IllegalArgumentException("Error: attempt to set orientation of plane with respect to an ill-defined vector");
        a = n.getX(0);
        b = n.getX(1);
        c = n.getX(2);
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
    public void setToCenter(Vector v) {
        v.setX(0,-d*a);
        v.setX(1,-d*b);
        v.setX(2,-d*c);
    }
    
    /**
     * Determines the point in the plane closest to the point x0 (first argument).
     * Uses second vector to hold result that is returned; if second vector is 
     * null, creates a new Vector instance.
     */
    public void setToNearestPoint(Vector x0, Vector point) {
        double factor = distanceTo(x0);
        point.setX(0,x0.getX(0)-factor*a);
        point.setX(1,x0.getX(1)-factor*b);
        point.setX(2,x0.getX(2)-factor*c);
    }

    /**
     * Perpendicular distance from the plane to the given point.
     */
    public double distanceTo(Vector x0) {
        return a*x0.getX(0) + b*x0.getX(1) + c*x0.getX(2) + d;
    }
    
    /**
     * Shifts the plane at fixed orientation so that it contains the given point.
     */
    public void moveTo(Vector r) {
        d = -(a*r.getX(0) + b*r.getX(1) + c*r.getX(2));
    }        
    
    /**
     * Dihedral angle between this plane and the given plane.
     */
    public double dihedralAngle(Plane p) {
        Vector n = space.makeVector();
        p.setToNormalVector(n);
        return a*n.getX(0) + b*n.getX(1) + c*n.getX(2);
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
    public boolean isPositiveSide(Vector p) {
        return distanceTo(p) > epsilon;
    }
    
    /**
     * Returns true if the given point is within a distance epsilon of this plane.
     */
    public boolean inPlane(Vector p) {
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
     * the given array of vectors.  Square is aligned so that its vertices
     * fall on the lines defined by the inPlaneVectors result.
     */
    public void setToInPlaneSquare(Vector x0, double size, Vector[] s) {
        setToInPlaneVectors(inPlane);
        Vector p1 = inPlane[0];
        Vector p2 = inPlane[1];
      //  System.out.println("work: "+x0.toString());
      //  System.out.println("p1 : "+p1.toString());
      //  System.out.println("p2 : "+p2.toString());
        double t = size/Math.sqrt(2.0);
        s[0].E(x0); s[0].PEa1Tv1(+t, p1);
        s[1].E(x0); s[1].PEa1Tv1(-t, p1);
        s[2].E(x0); s[2].PEa1Tv1(+t, p2);
        s[3].E(x0); s[3].PEa1Tv1(-t, p2);
    }
    
    /**
     * inPlaneSquare with square centered on point in plane closest to origin.
     */
    public void setToInPlaneSquare(double size, Vector[] s) {
        setToCenter(work0);
        setToInPlaneSquare(work0, size, s);
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
        Space space = Space3D.getInstance();
        Plane plane = new Plane(space, 2.8, 7.5, -2.1, 293.8);
    //    Plane plane = new Plane(5.0, 0.0, 0.0, 293.8);
        Vector[] p = new Vector[]{space.makeVector(),space.makeVector()};
        plane.setToInPlaneVectors(p);
        Vector n = space.makeVector();
        plane.setToNormalVector(n);
        System.out.println("Normal: "+n.toString());
        System.out.println("plane1: "+p[0].toString());
        System.out.println("plane2: "+p[1].toString());
        System.out.println("p1.n: "+p[0].dot(n));
        System.out.println("p2.n: "+p[1].dot(n));
        System.out.println("p1.p2: "+p[0].dot(p[1]));
    }
    
}//end of Plane
