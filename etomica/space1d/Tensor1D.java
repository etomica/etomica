package etomica.space1d;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Tensor1D implements etomica.space.Tensor {
    double xx;
    public static final Tensor1D ZERO = new Tensor1D();
    public static final Tensor1D IDENTITY = new Tensor1D(new double[] {1.0});
    public static final Tensor1D WORK = new Tensor1D();  //anything using WORK is not thread-safe
    public Tensor1D () {xx = 0.0;}
    public Tensor1D (double[] d) {this.E(d);}

    public int length() {return 1;}
    public double component(int i, int j) {return xx;}
    public void setComponent(int i, int j, double d) {xx = d;}
    public void E(Tensor1D t) {xx=t.xx;}
    public void E(Vector1D u1, Vector1D u2) {xx=u1.x*u2.x;}
    public void E(double a) {xx = a;}
    public void PE(Tensor1D t) {xx+=t.xx;}
    public void PE(int i, int j, double a) {xx += a;}
    public void PE(Vector1D u1, Vector1D u2) {xx+=u1.x*u2.x;}
    public double trace() {return xx;}
    
    public void E(etomica.space.Tensor t) {E((Tensor1D)t);}
    public void E(etomica.space.Vector u1, etomica.space.Vector u2) {E((Vector1D)u1, (Vector1D)u2);}
    public void PE(etomica.space.Tensor t) {PE((Tensor1D)t);}
    public void PE(etomica.space.Vector u1, etomica.space.Vector u2) {PE((Vector1D)u1, (Vector1D)u2);}
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
