/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space.Tensor;

/**
 * Implementation of a tensor for a 2-dimensional space.  Tensor is formed from four elements.
 */
public class Tensor2D implements etomica.space.Tensor, java.io.Serializable {

    double xx, xy, yx, yy;
    private static final long serialVersionUID = 1;
    
    public int D() {return 2;}

    /**
     * Default constructor sets all tensor elements to zero.
     */
    public Tensor2D () {xx = xy = yx = yy = 0.0;}
    
    /**
     * Constructs tensor with elements set by the given array.  Elements
     * are interpreted in order as xx, xy, yx, yy.
     */
    public Tensor2D (double[][] d) {
        this.E(d);
    }
    
    /**
     * Returns tensor in an array, writing rowwise, i.e., in order xx, xy, yx, yy.
     */
    public double[] toArray() {
        return new double[] {xx, xy, yx, yy};
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

    public double component(int i, int j) {
        return (i==0) ? ( (j==0) ? xx : xy ) : ( (j==0) ? yx : yy );
    }
    
    public void setComponent(int i, int j, double d) {
        if(i==0) {if(j==0) xx=d; else xy=d;}
        else     {if(j==0) yx=d; else yy=d;}
    }
    
    public void E(Tensor t) {
        xx=((Tensor2D)t).xx;
        xy=((Tensor2D)t).xy;
        yx=((Tensor2D)t).yx;
        yy=((Tensor2D)t).yy;
    }
    
    public void E(Vector[] v) {
        if(v.length != 2) {
            throw new IllegalArgumentException("Tensor requires 2 vectors to set its values");
        }
        xx = ((Vector2D)v[0]).x; xy = ((Vector2D)v[1]).x;
        yx = ((Vector2D)v[0]).y; yy = ((Vector2D)v[1]).y;
    }

    public void diagE(Vector v) {
        this.E(0.0);
        Vector2D v3 = (Vector2D)v;
        xx = v3.x;
        yy = v3.y;
    }
    
    public void assignTo(Vector[] v) {
        if(v.length != 2) {
            throw new IllegalArgumentException("Tensor requires 2 vectors for assignment");
        }
        ((Vector2D)v[0]).x = xx; ((Vector2D)v[1]).x = xy;
        ((Vector2D)v[0]).y = yx; ((Vector2D)v[1]).y = yy;
    }
   
    public void Ev1v2(Vector u1, Vector u2) {
        xx=((Vector2D)u1).x*((Vector2D)u2).x;
        xy=((Vector2D)u1).x*((Vector2D)u2).y;
        yx=((Vector2D)u1).y*((Vector2D)u2).x;
        yy=((Vector2D)u1).y*((Vector2D)u2).y;
    }
    
    public void E(double a) {
        xx = xy = yx = yy = a;
    }
    
    public void PE(double a) {
        xx+=a;
        xy+=a;
        yx+=a;
        yy+=a;
    }
    
    public void PE(Tensor t) {
        xx+=((Tensor2D)t).xx;
        xy+=((Tensor2D)t).xy;
        yx+=((Tensor2D)t).yx;
        yy+=((Tensor2D)t).yy;
    }
    
    public void PE(int i, int j, double a) {
        if(i==0) {if(j==0) xx+=a; else xy+=a;}
        else     {if(j==0) yx+=a; else yy+=a;}
    }
    
    public void PEv1v2(Vector v1, Vector v2) {
        xx+=((Vector2D)v1).x*((Vector2D)v2).x;
        xy+=((Vector2D)v1).x*((Vector2D)v2).y;
        yx+=((Vector2D)v1).y*((Vector2D)v2).x;
        yy+=((Vector2D)v1).y*((Vector2D)v2).y;
    }
    
    public void PEa1Tt1(double a1, Tensor t1) {
        xx += a1*((Tensor2D)t1).xx;
        xy += a1*((Tensor2D)t1).xy;
        yx += a1*((Tensor2D)t1).yx;
        yy += a1*((Tensor2D)t1).yy;
    }

    public void MEv1v2(Vector v1, Vector v2) {
        xx-=((Vector2D)v1).x*((Vector2D)v2).x;
        xy-=((Vector2D)v1).x*((Vector2D)v2).y;
        yx-=((Vector2D)v1).y*((Vector2D)v2).x;
        yy-=((Vector2D)v1).y*((Vector2D)v2).y;
    }
    
    public void ME(Tensor t) {
        xx-=((Tensor2D)t).xx;
        xy-=((Tensor2D)t).xy;
        yx-=((Tensor2D)t).yx;
        yy-=((Tensor2D)t).yy;
    }
    
    public double trace() {
        return xx + yy;
    }
    
    public void transpose(){
        double temp = xy; xy = yx; yx = temp;
    }
    
    public double determinant() {
        return xx*yy - xy*yx;
    }

    public void invert() {
        double det = determinant();
        double temp = xx; 
        xx = yy/det; yy = temp/det;
        xy = -xy/det; yx = -yx/det;
    }
    
    public void TE(double a) {
        xx*=a;  xy*=a; 
        yx*=a;  yy*=a;
    }
    
    public void TE(Tensor t){ Tensor2D u = (Tensor2D)t;
	    double txx=xx;double txy=xy;
	    double tyx=yx;double tyy=yy;
	    xx= txx*u.xx+txy*u.yx;   
	    xy= txx*u.xy+txy*u.yy; 
	    yx= tyx*u.xx+tyy*u.yx;   
	    yy= tyx*u.xy+tyy*u.yy; 
    }
    
    public void DE(Tensor t) {
        xx/=((Tensor2D)t).xx;
        xy/=((Tensor2D)t).xy;
        yx/=((Tensor2D)t).yx;
        yy/=((Tensor2D)t).yy;
    }
    
    public void E(double[] d) {
        if(d.length != 4) throw new IllegalArgumentException("Array size incorrect for tensor; (required, given): ("+4+", "+d);
        xx = d[0]; xy = d[1];
        yx = d[2]; yy = d[3];
    }
    
    public void E(double[][] d) {
        if(d.length != 2 || d[0].length != 2 || d[1].length != 2) throw new IllegalArgumentException("Array size incorrect for tensor; required: (2,2,2), given: ("+d.length+","+d[0].length+","+d[1].length+")");
        xx = d[0][0]; xy = d[0][1];
        yx = d[1][0]; yy = d[1][1];
    }
    
    public void assignTo(double[] d) {
        if(d.length != 4) throw new IllegalArgumentException("Array size incorrect for tensor; (required, given): ("+4+", "+d);
        d[0] = xx; d[1] = xy; 
        d[2] = yx; d[3] = yy;
    }
    
    public void assignTo(double[][] d) {
        d[0][0] = xx; d[0][1] = xy;
        d[1][0] = yx; d[1][1] = yy;
    }
    
    public boolean isNaN() {
        return Double.isNaN(xx) || Double.isNaN(xy) || Double.isNaN(yx) || Double.isNaN(yy);
    }
    
    public void map(IFunction f) {
        xx = f.f(xx);
        xy = f.f(xy);
        yx = f.f(yx);
        yy = f.f(yy);
    }
    
    public void transform(Vector v) {
        double x = xx * v.getX(0) + xy * v.getX(1);
        v.setX(1, yx * v.getX(0) + yy * v.getX(1));
        v.setX(0, x);
    }

    public boolean equals(Tensor t) {
        Tensor2D A = (Tensor2D)t;
        return A.xx==xx && A.xy==xy && A.yx==yx && A.yy==yy;
    }

   public String toString() {
        return "("+xx+", "+xy+")\n("+yx+", "+yy+")";
    }
}
