package etomica.space;

import etomica.NearestImageTransformer;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.utility.Function;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public abstract class Vector implements java.io.Serializable, Cloneable {

    /**
     * Support of implementation of Cloneable interface. Returns a new Vector
     * with elements equal to this one.
     */
    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException ex) {
            throw new InternalError(ex.toString());
        }
    }

    /**
     * Dimension of the space occupied by the vector. Number of element in the
     * vector.
     */
    public abstract int D();

    /**
     * Assigns the components of this vector to the elements of the given array.
     * Does not check that array length is consistent with vector dimension.
     * Inverse of the E method.
     */
    public abstract void assignTo(double[] array);

    /**
     * Returns a new array with elements equal to the components of this vector.
     */
    public abstract double[] toArray();

    /**
     * Returns true if all corresponding elements of this and the given vector
     * are equal; returns false otherwise.
     */
    public abstract boolean equals(Vector v);

    /**
     * Computes the spherical-coordinate representation of this vector, and
     * returns them in the given array.
     */
    public abstract void sphericalCoordinates(double[] result);

    /**
     * Computes the spherical-coordinate representation of this vector, and
     * returns them in the given array.
     */
    public double[] toSphericalCoordinateArray() {
        double[] array = new double[D()];
        sphericalCoordinates(array);
        return array;
    }

    /**
     * Vector components corresponding to the given index. For example,
     * x-component is given for i = 0. If index exceeds dimension of vector, no
     * error is given and the last element of the vector is returned.
     */
    public abstract double x(int i);

    /**
     * Sets the i-th component of the vector to the given value d. If index
     * exceeds dimension of vector, no error is given and the last element of
     * the vector is assigned.
     */
    public abstract void setX(int i, double d); 
    
    /**
     * Sets the components of this vector equal to those of the given vector.
     */
    public abstract void E(Vector u);

    /**
     * Sets all components of this vector equal to the given value.
     */
    public abstract void E(double a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given array. Inverse of the assignTo method.
     */
    public abstract void E(double[] a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given integer array.
     */
    public abstract void E(int[] a);

    /**
     * Plus-equals (+=) operation. Increments each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void PE(Vector u);

    /**
     * Plus-equals (+=) operation applied to one component.
     * 
     * @param i
     *            index of the component being incremented
     * @param a
     *            increment value
     */
    public abstract void PE(int i, double a);

    /**
     * Plus-equals (+=) operation, causing a constant value to be added to all
     * components of this vector.
     */
    public abstract void PE(double a); //adds a constant value to all elements

    /**
     * Minus-equals (-=) operation. Decrements each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void ME(Vector u);

    /**
     * Times-equals (*=) operation. Multiplies each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void TE(Vector u);

    /**
     * Divide-equals (/=) operation. Divides each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void DE(Vector u);

    /**
     * Multiplies all components by a constant.
     */
    public abstract void TE(double a); //multipies all components by the

    // constant a

    /**
     * Multiplies one component, of index i, by a constant.
     */
    public abstract void TE(int i, double a);

    /**
     * Operation (= a1 * v1); sets the components of this to those of the given
     * vector multiplied by the given constant.
     */
    public abstract void Ea1Tv1(double a, Vector v1); //sets this vector to a*u

    /**
     * Sets this vector equal to (v1 + a1 * v2).
     */
    public abstract void Ev1Pa1Tv2(Vector v1, double a1, Vector v2);

    /**
     * Increments (+=) components by (a1 * v1)
     */
    public abstract void PEa1Tv1(double a, Vector v1);

    /**
     * Sets the components of this vector equal to (v1 + v2)
     */
    public abstract void Ev1Pv2(Vector v1, Vector v2); //sets equal to sum of

    /**
     * Adds/subtracts a scalar to each element, with the choice of 
     * addition or subtraction for each element determined by the corresponding
     * element of the vector v1.  In symbols, this is += a1*sgn(v1).  Zero
     * is added if the corresponding element of v1 is zero.
     */
    public abstract void PEa1SGNv1(double a1, Vector v1);
    
    // vectors
    /**
     * Sets the components of this vector equal to (v1 - v2)
     */
    public abstract void Ev1Mv2(Vector v1, Vector v2);

    /**
     * Adds (Plus, +) this vector to the given one and returns sum in a new vector.
     */
    public abstract Vector P(Vector u);

    /**
     * Subtracts (Minus, -) given vector from this returns result in a new vector.
     */
    public abstract Vector M(Vector u);

    /**
     * Multiplies (Times, *) given vector times this and returns result 
     * in a new vector.
     */
    public abstract Vector T(Vector u); 

    /**
     * Divides (/) this vector by the given one and returns result in a new vector.
     */
    public abstract Vector D(Vector u);

    /**
     * Returns square of vector resulting from subtracting given vector
     * from this one.  Neither vector is changed by this operation.
     */
    public abstract double Mv1Squared(Vector v1);

    /**
     * Replaces each component of this vector with its absolute value.
     */
    public abstract void abs(); //replaces each component with its absolute

    /**
     * Replaces each component of this vector with its value modulo a
     */
    public abstract void mod(double a); 

    /**
     * Replaces each component of this vector with its value modulo
     * the corresponding component in u
     */
    public abstract void mod(Vector u); //each component replaced with itself

    /**
     * Sets this equal to (r mod u) - r
     */
    public abstract void EModShift(Vector r, Vector u);

    /**
     * Hard to explain.
     */
    public abstract void EMod2Shift(Vector r, Vector u);

    /**
     * Returns the minimum of all components of this vector.
     */
    public abstract double min();

    /**
     * Returns the maximum of all components of this vector.
     */
    public abstract double max();

    /**
     * Returns the square magnitude of this vector, e.g., x^2 + y^2 for D = 2.
     */
    public abstract double squared();

    /**
     * Normalizes this vector, so this.squared() == 1.  Divides all 
     * components by Math.sqrt(this.squared()).
     */
    public abstract void normalize();

    /**
     * Returns true if all components of this vector are zero; false otherwise.
     */
    public abstract boolean isZero();

    /**
     * Returns the dot product of this vector with the given one.
     */
    public abstract double dot(Vector u); 

    /**
     * Returns the 3D cross product of this vector with a 2D vector.
     */
    public abstract etomica.space3d.Vector3D cross(etomica.space2d.Vector2D u);

    /**
     * Cross product of this vector with the given vector
     */
    public abstract etomica.space3d.Vector3D cross(etomica.space3d.Vector3D u); 

    /**
     * Replaces this vector with its cross product with the given vector.
     */
    public abstract void XE(etomica.space3d.Vector3D u); //replaces this vector

    /**
     * Applies the given tensor transformation to this vector, replaced its
     * elements with the transformed values.
     */
    public abstract void transform(Tensor A);

    /**
     * Assigns each components to (its own) random value between 0 and d.
     */
    public abstract void setRandom(double d);

    /**
     * Assigns each component to (its own) random value between 0 and v(i)
     */
    public abstract void setRandom(Vector v);
    
    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public abstract void setRandomSphere();


    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public abstract void setRandomCube();

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public abstract void setRandomInSphere();

    /**
     * Assigns this vector to a new one chosen randomly in the cube of
     * size 2*d (selecting +/- d) centered on its present position.
     */
    public abstract void randomStep(double d); 
    
    /**
     * Adds several terms of the form a*u to this vector.  Adds the sum
     * of a[i]*u[i] to the i-th component component, doing this to all
     * components.  Does not check compatibility of array lengths.
     */
    public final void PEa1Tv1(double[] a, Vector[] u) { 
        for (int i = a.length - 1; i >= 0; i--) {
            PEa1Tv1(a[i], u[i]);
        }
    }

    /**
     * Returns a new vector equal to cross product of this vector 
     * with the given one
     */
    public etomica.space3d.Vector3D cross(Vector u) {
        if (u instanceof Vector3D) {
            return cross((Vector3D) u);
        } else if (u instanceof Vector2D) {
            return cross((Vector2D) u);
        } else
            return null;
    }

    /**
     * Performs a random rotation of this vector within a cone of 
     * angle theta on its current position.
     */
    public abstract void randomRotate(double thetaStep);
    
    /**
     * Returns true if any element of the vector is not-a-number.
     */
    public abstract boolean isNaN();
    
    /**
     * Applies the given function to each element of the vector.
     */
    public abstract void map(Function f);

    /**
     * Applies the given tensor transformation to the nearest-image transform
     * of this minus r0.  Specifically, performs<code><p>
     * this.ME(r0);<br>
     * nit.nearestImage(this);<br>
     * this.transform(A);<br>
     * this.PE(r0)</code>
     */
    public void transform(NearestImageTransformer nit, Vector r0, Tensor A) {
        this.ME(r0);
        nit.nearestImage(this);
        this.transform(A);
        this.PE(r0);
    }

}