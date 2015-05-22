/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;

/**
 * Quaternion object , its properties and associated operations
 *
 * @author rsubrama
 */
// TODO : Test all methods
public class Quaternion extends Object {
	protected double q0,q1,q2,q3;
	protected static final ISpace space = Space3D.getInstance();
	protected final IVector vec = space.makeVector();
	protected double norm;	
	
	/**
       Constructs the Quaternion q = p[0] + p[1]*i + p[2]*j + p[3]*k
       @author rsubrama
       @return Quaternion q
	 */
	public Quaternion(double [] p) {
	    if (p.length != 4) throw new RuntimeException("Length of the double array should be 4");
	    q0 = p[0];
	    q1 = p[1];
	    q2 = p[2];
	    q3 = p[3];	    
	    updateQuaternion();	    
	}
	
	/**
       Constructs the Quaternion q = a + b
       @author rsubrama
       @return Quaternion q
	 */
	public Quaternion(double a, IVectorMutable b) {
        q0 = a;
        q1 = b.getX(0);
        q2 = b.getX(1);
        q3 = b.getX(2);
        updateQuaternion();
    }
	
	/**
	 * Constructs the unit Quaternion q = cos(angle/2) + axis*sin(angle/2)
	 * @author rsubrama
	 * @param axis
	 * @param angle
	 * @return Quaternion q
	 */
	public Quaternion(IVectorMutable axis, double angle) {
	    axis.normalize();
        double halfAngle = angle/2.0;
        q0 = Math.cos(halfAngle);
        q1 = axis.getX(0)*Math.sin(halfAngle);
        q2 = axis.getX(1)*Math.sin(halfAngle);
        q3 = axis.getX(2)*Math.sin(halfAngle);
        updateQuaternion();
	}
	
	/**
       Constructs the Quaternion q = a + b*i + c*j + d*k
       @author rsubrama
       @return Quaternion q
	 */
	public Quaternion(double a, double b, double c, double d) {
	    q0 = a;
	    q1 = b;
	    q2 = c;
	    q3 = d;
	    updateQuaternion();
	}
	
	/**
       Constructs the Quaternion q = a + b[0]*i + b[1]*j + b[2]*k
       @author rsubrama
       @return Quaternion q
	 */
	public Quaternion(double a, double [] b) {
	    q0 = a;
	    if (b.length != 3) throw new RuntimeException("Length of double array should be 3");
	    q1 = b[0];
	    q2 = b[1];
	    q3 = b[2];
	    updateQuaternion();
	}
	
	/**
	 * Constructs the Quaternion q for an input OrientationFull3D with two vectors.
	 * First creates the three basis vectors from the orientation. Then converts that 
	 * into a rotation matrix and finally into a Quaternion.
	 * This Quaternion represents the rotation of the {i,j,k} coordinate frame to the 
     * orthonormal coordinate frame formed by the the basis vectors.
	 * @author rsubrama
	 * @param or
	 * @return Quaternion q
	 */
	public Quaternion(OrientationFull3D or) {
        // Get basis vectors from orientation
        IVectorMutable [] v = space.makeVectorArray(3);
        IVectorMutable ex = space.makeVector();
        IVectorMutable ey = space.makeVector();
        IVectorMutable ez = space.makeVector();
        ex.E(0);
        ey.E(0);
        ez.E(0);
        ex.setX(0, 1);
        ey.setX(1, 1);
        ez.setX(2, 1);
        v[0].E(or.getDirection());
        v[1].E(or.getSecondaryDirection());
        v[2].E(v[0]);
        v[2].XE(v[1]);
        v[2].normalize();
        // Direction cosines to Rotation matrix
        double a00 = v[0].dot(ex);
        double a01 = v[0].dot(ey);
        double a02 = v[0].dot(ez);
        double a10 = v[1].dot(ex);
        double a11 = v[1].dot(ey);
        double a12 = v[1].dot(ez);
        double a20 = v[2].dot(ex);
        double a21 = v[2].dot(ey);
        double a22 = v[2].dot(ez);
        double tr = a00 + a11 + a22;
        
        // Rotation matrix to quaternion. Code taken from:
        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        if (tr > 0) {
            double S = Math.sqrt(tr+1.0) * 2; // S=4*q0
            q0 = 0.25 * S;
            q1 = (a21 - a12) / S;
            q2 = (a02 - a20) / S;
            q3 = (a10 - a01) / S;
        } else if ((a00 > a11)&&(a00 > a22)) {
            double S = Math.sqrt(1.0 + a00 - a11 - a22) * 2; // S=4*q1
            q0 = (a21 - a12) / S;
            q1 = 0.25 * S;
            q2 = (a01 + a10) / S;
            q3 = (a02 + a20) / S;
        } else if (a11 > a22) {
            double S = Math.sqrt(1.0 + a11 - a00 - a22) * 2; // S=4*q2
            q0 = (a02 - a20) / S;
            q1 = (a01 + a10) / S;
            q2 = 0.25 * S;
            q3 = (a12 + a21) / S;
        } else {
            double S = Math.sqrt(1.0 + a22 - a00 - a11) * 2; // S=4*q3
            q0 = (a10 - a01) / S;
            q1 = (a02 + a20) / S;
            q2 = (a12 + a21) / S;
            q3 = 0.25 * S;
        }
        updateQuaternion();
    }
	
	/**
	 * Constructs the Quaternion q for an input basis vector array. Converts the basis vector array 
     * into a rotation matrix and finally into a Quaternion.
     * This Quaternion represents the rotation of the {i,j,k} coordinate frame to the
     * orthonormal coordinate frame formed by the the basis vectors.
     * @author rsubrama
	 * @param basisVec
	 */
	public Quaternion(IVector[] basisVec) {        
        IVectorMutable ex = space.makeVector();
        IVectorMutable ey = space.makeVector();
        IVectorMutable ez = space.makeVector();
        ex.E(0);
        ey.E(0);
        ez.E(0);
        ex.setX(0, 1);
        ey.setX(1, 1);
        ez.setX(2, 1);
        
        // Direction cosines to Rotation matrix 
        double a00 = basisVec[0].dot(ex);
        double a01 = basisVec[0].dot(ey);
        double a02 = basisVec[0].dot(ez);
        double a10 = basisVec[1].dot(ex);
        double a11 = basisVec[1].dot(ey);
        double a12 = basisVec[1].dot(ez);
        double a20 = basisVec[2].dot(ex);
        double a21 = basisVec[2].dot(ey);
        double a22 = basisVec[2].dot(ez);
        double tr = a00 + a11 + a22;
        
        // Rotation matrix to quaternion. Code taken from: 
        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        if (tr > 0) {
            double S = Math.sqrt(tr+1.0) * 2; // S=4*q0 
            q0 = 0.25 * S;
            q1 = (a21 - a12) / S;
            q2 = (a02 - a20) / S;
            q3 = (a10 - a01) / S;
        } else if ((a00 > a11) && (a00 > a22)) {
            double S = Math.sqrt(1.0 + a00 - a11 - a22) * 2; // S=4*q1 
            q0 = (a21 - a12) / S;
            q1 = 0.25 * S;
            q2 = (a01 + a10) / S;
            q3 = (a02 + a20) / S; 
        } else if (a11 > a22) {
            double S = Math.sqrt(1.0 + a11 - a00 - a22) * 2; // S=4*q2
            q0 = (a02 - a20) / S;
            q1 = (a01 + a10) / S;
            q2 = 0.25 * S;
            q3 = (a12 + a21) / S;
        } else {
            double S = Math.sqrt(1.0 + a22 - a00 - a11) * 2; // S=4*q3
            q0 = (a10 - a01) / S;
            q1 = (a02 + a20) / S;
            q2 = (a12 + a21) / S;
            q3 = 0.25 * S;
        }
        updateQuaternion();
    }
	
	/**
	 * Constructs the Quaternion q given an input rotation matrix. Note that the determinant of the 
	 * input matrix should be +1. Otherwise it doesn't represent a rotation
	 * @author rsubrama
	 * @param rotMatrix
	 * @return Quaternion q
	 */
	public Quaternion(Jama.Matrix rotMatrix) {
	    if (rotMatrix.det() != 1 && (rotMatrix.det() - 1)*(rotMatrix.det() - 1) > 1E-15) throw new RuntimeException("Determinant of matrix not equal to +1.\nTherefore it doesn't represent a rotation");
	    double [][] a = rotMatrix.getArray();
        double tr = a[0][0] + a[1][1] + a[2][2];
        
        // Rotation matrix to quaternion. Code taken from: 
        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        if (tr > 0) {
            double S = Math.sqrt(tr+1.0) * 2; // S=4*q0 
            q0 = 0.25 * S;
            q1 = (a[2][1] - a[1][2]) / S;
            q2 = (a[0][2] - a[2][0]) / S;
            q3 = (a[1][0] - a[0][1]) / S;
        } else if ((a[0][0] > a[1][1])&&(a[0][0] > a[2][2])) {
            double S = Math.sqrt(1.0 + a[0][0] - a[1][1] - a[2][2]) * 2; // S=4*q1 
            q0 = (a[2][1] - a[1][2]) / S;
            q1 = 0.25 * S;
            q2 = (a[0][1] + a[1][0]) / S;
            q3 = (a[0][2] + a[2][0]) / S; 
        } else if (a[1][1] > a[2][2]) {
            double S = Math.sqrt(1.0 + a[1][1] - a[0][0] - a[2][2]) * 2; // S=4*q2
            q0 = (a[0][2] - a[2][0]) / S;
            q1 = (a[0][1] + a[1][0]) / S;
            q2 = 0.25 * S;
            q3 = (a[1][2] + a[2][1]) / S;
        } else {
            double S = Math.sqrt(1.0 + a[2][2] - a[0][0] - a[1][1]) * 2; // S=4*q3
            q0 = (a[1][0] - a[0][1]) / S;
            q1 = (a[0][2] + a[2][0]) / S;
            q2 = (a[1][2] + a[2][1]) / S;
            q3 = 0.25 * S;
        }
        updateQuaternion();
    }
	
    /**
     * First checks if all the components are numbers. Then recomputes the norm and the vector
     * of the Quaternion that has been operated upon.
     * @author rsubrama
     */
    protected void updateQuaternion() {
        double [] q = this.getDoubleArray();
        for (int i = 0; i<q.length; i++) {
            if (Double.isInfinite(q[i]) || Double.isNaN(q[i])) throw new RuntimeException("Oops! q["+i+"] = "+q[i]);
        }
        ((IVectorMutable) vec).setX(0, q1);
        ((IVectorMutable) vec).setX(1, q2);
        ((IVectorMutable) vec).setX(2, q3);
        norm = Math.sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    }	
	
	// properties
	/**
	 * Returns the q0 component of the Quaternion
	 * @author rsubrama
	 * @return double q0
	 */
	public double getQ0() {
	    return q0;
	}
	
	/**
	 * Returns the q1 component of the Quaternion
	 * @author rsubrama
	 * @return double q1
	 */
	public double getQ1() {
        return q1;
    }
	
	/**
	 * Returns the q2 component of the Quaternion
	 * @author rsubrama
	 * @return double q2
	 */
	public double getQ2() {
        return q2;
    }
	
	/**
	 * Returns the q3 component of the Quaternion
	 * @author rsubrama
	 * @return double q3
	 */
	public double getQ3() {
        return q3;
    }
	
	/**
	 *  Returns the scalar part of the Quaternion
	 *  @author rsubrama
	 *  @return double q0
	 */
	public double getScalar() {
	    return q0;
	}
	
	/**
	 * Returns the vector part of the Quaternion as IVectorMutable
	 * @author rsubrama
	 * @return IVectorMutable v = q1*i + q2*j + q3*k
	 */
	public IVectorMutable getVectorMutable() {
	    return (IVectorMutable) vec;
	}
	
	/**
	 * Returns the vector part of the Quaternion as IVector
	 * @author rsubrama
	 * @return IVector v = q1*i + q2*j + q3*k
	 */
	public IVector getVector() {
	    return vec;
	}
	
	/**
	 * Assigns the angle and the axis associated with the Quaternion. Note, axis is a unit vector
	 * @author rsubrama
	 * @param axis 
	 * @param angle 
	 */
	public void getAxisAngle(IVectorMutable axis, double angle) {
	    if (! this.isUnitQuaternion()) throw new RuntimeException("This is not a unit quaternion and therefore not a rotation operator. Normalize it first");
	    if (q0 > 1.0) q0 = 1.0;
	    if (q0 < -1.0) q0 = -1.0;
	    angle = 2*Math.acos(q0);
	    axis.E(vec);
	    axis.normalize();
	}
	
	/**
	 * Returns the vector associated with the Quaternion as a double array
	 * @author rsubrama
	 * @return double [] v
	 */
	public double[] getVectorAsDoubleArray() {
	    double [] v = new double[3];
	    vec.assignTo(v);
	    return v;
	}
	
	/**
	 * Returns the magnitude of the vector associated with the Quaternion
	 * @author rsubrama
	 * @return double magnitude = |vec|
	 */
	public double getVectorMagnitude() {
	    return Math.sqrt(vec.squared());
	}
	
	/**
	 * Returns the Quaternion as a double array
	 * @author rsubrama
	 * @return double [] q
	 */
	public double[] getDoubleArray() {
	    double [] q = {q0,q1,q2,q3};
	    return q;
	}
	
	/**
	 * Returns the rotation matrix associated with the Quaternion q as a Jama.Matrix
	 * @author rsubrama
	 * @return Jama.Matrix A
	 */
	public Jama.Matrix getRotationMatrix() {
	    double[][] a = new double[3][3];
	    a[0][0] = 2*q0*q0 - 1 + 2*q1*q1;
	    a[0][1] = 2*q1*q2 - 2*q0*q3;
	    a[0][2] = 2*q1*q3 + 2*q0*q2;
	    
	    a[1][0] = 2*q1*q2 + 2*q0*q3;
	    a[1][1] = 2*q0*q0 - 1 + 2*q2*q2;
	    a[1][2] = 2*q2*q3 - 2*q0*q1;
	    
	    a[2][0] = 2*q1*q3 - 2*q0*q2;
	    a[2][1] = 2*q2*q3 + 2*q0*q1;
	    a[2][2] = 2*q0*q0 - 1 + 2*q3*q3;
	    Jama.Matrix A = new Jama.Matrix(a);
	    return A;
	}
	
	/**
	 * Returns an orientation corresponding to the rotation of cartesian coordinate frame {i,j,k}
	 * by the axis and angle associated with the Quaternion q.
	 * @author rsubrama
	 * @return OrientationFull3D or
	 */
	public OrientationFull3D getOrientation() {
	    OrientationFull3D or = (OrientationFull3D) space.makeOrientation();
	    IVectorMutable [] e = space.makeVectorArray(3);
	    for (int i=0; i<3; i++) {
	        e[i].E(0);
	        e[i].setX(i, 1);
	    }
	    or.setDirections(e[0], e[1]);
	    this.rotateOrientation(or);
	    return or;
	}
	
	/**
	 * Returns whether the Quaternion is pure or not
	 * @author rsubrama
	 * @return boolean true/false
	 */
	public boolean isPureQuaternion() {
	    if (q0 == 0 || q0*q0 < 1E-15) return true;
	    return false;
	}
	
	/**
	 * Returns whether the Quaternion is equal to another given Quaternion p
	 * @author rsubrama
	 * @param p
	 * @return boolean true/false
	 */
	public boolean isEqualTo(Quaternion p) {
	    if (q0 == p.q0 && q1 == p.q1 && q2 == p.q2 && q3 == p.q3) return true;
	    if ((q0 - p.q0)*(q0 - p.q0) < 1E-15 && (q1 - p.q1)*(q1 - p.q1) < 1E-15 && (q2 - p.q2)*(q2 - p.q2) < 1E-15 && (q3 - p.q3)*(q3 - p.q3) < 1E-15) return true;
	    return false;
	}
	
	/**
	 * Returns whether the norm of the Quaternion is 1 or not
	 * @author rsubrama
	 * @return boolean true/false
	 */
	public boolean isUnitQuaternion() {
	    if (norm == 1 || (norm - 1)*(norm - 1) < 1E-15) return true;
	    return false;
	}
	
	// operations
	/**
	 * Adds the input Quaternion components to its components. Addition is done component wise 
	 * @author rsubrama
	 * @param p
	 */
	public void PE(Quaternion p) {
	    q0 += p.q0;
	    q1 += p.q1;
	    q2 += p.q2;
	    q3 += p.q3;
	    updateQuaternion();
	}
	
	/**
	 * Multiplies all the components of the Quaternion q by the given scalar
	 * @author rsubrama
	 * @param a 
	 */
	public void TE(double a) {
	    if (q0 != 0) q0 *= a;
	    if (q1 != 0) q1 *= a;
	    if (q2 != 0) q2 *= a;
	    if (q3 != 0) q3 *= a;
	    updateQuaternion();
	}
	
	/**
	 * Returns the Quaternion r resulting from multiplying the input Quaternion p to the current Quaternion q.
	 * Can be thought of as a pre-multiplication with the current Quaternion.
	 * Mulitplication is done based on Hamilton's rules and result in (r = pq):
	 * r = p0*q0 - pVec.dot(qVec) + p0*qVec + q0*pVec + pVec X qVec
	 * @author rsubrama 
	 * @param p
	 * @return Quaternion r = pq
	 */
	public Quaternion multiply(Quaternion p) {
	    double a = p.q0*q0 - p.vec.dot(vec);
	    IVectorMutable newVec = space.makeVector();
	    newVec.E(p.vec);
	    newVec.XE(vec); 
	    newVec.PEa1Tv1(q0, p.vec);
	    newVec.PEa1Tv1(p.q0, vec);
	    return new Quaternion(a, newVec);
	}
	
	/**
	 * Returns the complex conjugate of the Quaternion q.
	 * @author rsubrama
	 * @return Quaternion q* = q0 - qVec
	 */
	public Quaternion conjugate() {
	    IVectorMutable vNew = space.makeVector();
	    vNew.E(vec);
	    vNew.TE(-1);
	    for (int i=0; i<3; i++) {
	        if (vNew.getX(i) == 0) vNew.setX(i, 0); 
	    }
	    return new Quaternion(q0, vNew);
	}
	
	/**
	 * Returns the norm of the Quaternion q defined as N = sqrt(q*q)
	 * @author rsubrama
	 * @return double N = sart(q0^2 + q1^2 + q2^2 + q3^2)
	 */
	public double getNorm(){
	    return norm;
	}
	
	/**
	 * Normalizes the Quaternion q i.e. divides all components by the norm
	 * @author rsubrama
	 */
	public void normalize(){
	    if (norm == 0) throw new RuntimeException("Norm is zero. Can't normalize this Quaternion");
	    q0 /= norm;
	    q1 /= norm;
	    q2 /= norm;
	    q3 /= norm;
	    updateQuaternion();
	}
	
	/**
	 * Returns the inverse of the Quaternion q given by q^-1 = 1/(norm^2) * q*
	 * @author rsubrama  
	 * @return Quaternion q^-1
	 */
	public Quaternion inverse(){
	    if (norm == 0) throw new RuntimeException("Norm is zero. Inverse doesn't exist");
	    Quaternion p = this.conjugate();
	    p.TE(1/(this.norm*this.norm));
	    return p;
	}
	
	/**
	 * Rotates the input vector by using the Quaternion rotation operator L(v) = qvq*.
	 * w = L(v) = qvq* = (2*q0^2 -1)v + 2*(qVec.dot(v))*qVec + 2*q0*(qVec X v)
	 * @author rsubrama
	 * @param v
	 */
	public void rotateVector(IVectorMutable v) {	    
	    IVectorMutable w = space.makeVector();
	    w.E(v);
	    v.E(vec);
	    v.XE(w);
	    v.TE(2*q0);
	    v.PEa1Tv1(2*vec.dot(w), vec);
	    v.PEa1Tv1((2*q0*q0 - 1), w);	    
	}
	
	/**
	 * Rotates the input orientation by using the Quaternion rotation operator L(v) = qvq*.
	 * Rotates the two vectors that form the orientation.
	 * @author rsubrama 
	 * @param or
	 */
	public void rotateOrientation(OrientationFull3D or) {
	    this.rotateVector((IVectorMutable)or.getDirection());
	    this.rotateVector((IVectorMutable)or.getSecondaryDirection());
	}
	
	/**
	 * Returns the string representation of the Quaternion q
	 * @author rsubrama
	 * @return String
	 */
	public String toString() {
	    String str = q0+" ";
	    double [] q = this.getDoubleArray();
	    for (int i=1; i<q.length; i++){
	        String str1 = "";
	        if (i == 1) str1 = "i";
	        if (i == 2) str1 = "j";
	        if (i == 3) str1 = "k";
	        if (q[i] >= 0) {
	            str += "+ "+q[i]+str1+" "; 
	        }
	        else {
	            str += "- "+(-q[i])+str1+" ";
	        }
	    }
	    return str;
	}
}
