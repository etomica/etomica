/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space1d;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space.Tensor;

/**
 * Implementation of a tensor for a 1-dimensional space. In this case the tensor
 * is a trivial object formed from just one element.
 */
public class Tensor1D implements etomica.space.Tensor, java.io.Serializable {

    double xx;
    private static final long serialVersionUID = 1;

    /**
     * Default constructor sets tensor element to zero.
     */
    public Tensor1D() {
        xx = 0.0;
    }

    public Tensor1D(double d) {
        this.E(d);
    }

    public double[] toArray() {
        return new double[] { xx };
    }

    /**
     * Support of implementation of Cloneable interface. Returns a new Tensor
     * with elements equal to this one.
     */
    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException ex) {
            throw new InternalError(ex.toString());
        }
    }

    public int D() {
        return 1;
    }

    public double component(int i, int j) {
        return xx;
    }

    public void setComponent(int i, int j, double d) {
        xx = d;
    }

    public void E(double a) {
        xx = a;
    }
    
    public void E(double[][] a) {
        xx = a[0][0];
    }
    
    public void E(Vector[] v) {
        if(v.length != 1) {
            throw new IllegalArgumentException("Tensor requires 1 vector to set its values");
        }
        xx = ((Vector1D)v[0]).x;
    }

    public void diagE(Vector v) {
        xx = ((Vector1D)v).x;
    }
    
    public void assignTo(Vector[] v) {
        if(v.length != 1) {
            throw new IllegalArgumentException("Tensor requires 1 vector to set its values");
        }
        ((Vector1D)v[0]).x = xx;
        
    }

    public void PE(double a) {
        xx += a;
    }

    public void PE(int i, int j, double a) {
        xx += a;
    }

    public double trace() {
        return xx;
    }

    public void transpose() {}

    public void invert() {
        xx = 1.0 / xx;
    }

    public void E(Tensor t) {
        xx = ((Tensor1D) t).xx;
    }

    public void Ev1v2(Vector u1, Vector u2) {
        xx = ((Vector1D) u1).x * ((Vector1D) u2).x;
    }

    public void PE(Tensor t) {
        xx += ((Tensor1D) t).xx;
    }

    public void PEv1v2(Vector u1, Vector u2) {
        xx += ((Vector1D) u1).x * ((Vector1D) u2).x;
    }
    
    public void PEa1Tt1(double a1, Tensor t1) {
        xx += a1*((Tensor1D)t1).xx;
    }

    public void MEv1v2(Vector u1, Vector u2) {
        xx -= ((Vector1D) u1).x * ((Vector1D) u2).x;
    }

    public void ME(Tensor t) {
        xx -= ((Tensor1D) t).xx;
    }

    public void TE(double a) {
        xx *= a;
    }

    public void TE(Tensor t) {
        xx *= ((Tensor1D) t).xx;
    }

    public void DE(Tensor t) {
        xx /= ((Tensor1D) t).xx;
    }

    public void E(double[] d) {
        if (d.length != 1)
            throw new IllegalArgumentException("Array size incorrect for tensor");
        xx = d[0];
    }

    public void assignTo(double[] d) {
        if (d.length != 1)
            throw new IllegalArgumentException("Array size incorrect for tensor");
        d[0] = xx;
    }

    public void assignTo(double[][] d) {
        d[0][0] = xx;
    }

    public boolean isNaN() {
        return Double.isNaN(xx);
    }

    public void map(IFunction f) {
        xx = f.f(xx);
    }
    
    public void transform(Vector v) {
        v.TE(xx);
    }

    public double determinant() {
        return xx;
    }

    public boolean equals(Tensor t) {
        return xx == ((Tensor1D)t).xx;
    }
}
