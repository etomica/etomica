package etomica.space1d;

import etomica.Space1D;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Tensor implements Tensor {
    double xx;
    public static final Tensor ZERO = new Tensor();
    public static final Tensor IDENTITY = new Tensor(new double[] {1.0});
    public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
    public Tensor () {xx = 0.0;}
    public Tensor (double[] d) {this.E(d);}

    public int length() {return Space1D.D;}
    public double component(int i, int j) {return xx;}
    public void setComponent(int i, int j, double d) {xx = d;}
    public void E(Tensor t) {xx=t.xx;}
    public void E(Vector u1, Vector u2) {xx=u1.x*u2.x;}
    public void E(double a) {xx = a;}
    public void PE(Tensor t) {xx+=t.xx;}
    public void PE(int i, int j, double a) {xx += a;}
    public void PE(Vector u1, Vector u2) {xx+=u1.x*u2.x;}
    public double trace() {return xx;}
    
    public void E(Tensor t) {E((Tensor)t);}
    public void E(Vector u1, Vector u2) {E((Vector)u1, (Vector)u2);}
    public void PE(Tensor t) {PE((Tensor)t);}
    public void PE(Vector u1, Vector u2) {PE((Vector)u1, (Vector)u2);}
    public void TE(double a) {xx*=a;}
    public void E(double[] d) {
        if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
        xx = d[0];
    }
    public void assignTo(double[] d) {
        if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
        d[0] = xx;
    }
    
}