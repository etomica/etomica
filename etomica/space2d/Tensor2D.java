package etomica.space2d;

import etomica.space.Tensor;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Tensor2D implements etomica.space.Tensor {
    public int length() {return 2;}
    double xx, xy, yx, yy;
    public static final Tensor2D ZERO = new Tensor2D();
    public static final Tensor2D WORK = new Tensor2D();  //anything using WORK is not thread-safe
    public Tensor2D () {xx = xy = yx = yy = 0.0;}
    public Tensor2D (double[] d) {
        this.E(d);
    }

    public double component(int i, int j) {
        return (i==0) ? ( (j==0) ? xx : xy ) : ( (j==0) ? yx : yy );}
    public void setComponent(int i, int j, double d) {
        if(i==0) {if(j==0) xx=d; else xy=d;}
        else     {if(j==0) yx=d; else yy=d;}
    }
    public void E(Tensor2D t) {xx=t.xx; xy=t.xy; yx=t.yx; yy=t.yy;}
    public void E(Vector2D u1, Vector2D u2) {xx=u1.x*u2.x; xy=u1.x*u2.y; yx=u1.y*u2.x; yy=u1.y*u2.y;}
    public void E(double a) {xx = xy = yx = yy = a;}
    public void PE(Tensor2D t) {xx+=t.xx; xy+=t.xy; yx+=t.yx; yy+=t.yy;}
    public void PE(int i, int j, double a) {
        if(i==0) {if(j==0) xx+=a; else xy+=a;}
        else     {if(j==0) yx+=a; else yy+=a;}
    }
    public void PE(Vector2D u1, Vector2D u2) {xx+=u1.x*u2.x; xy+=u1.x*u2.y; yx+=u1.y*u2.x; yy+=u1.y*u2.y;}
    public double trace() {return xx + yy;}
    public void transpose(){
    	double temp = 0.0;
    	temp = xy; xy = yx; yx = temp;
    }
    public void inverse() {
    	double det = xx*yy -xy*yx;
    	double temp =0.0;
    	temp = xx; xx = yy/det; yy = temp/det;
    	temp = xy; xy = -yx/det; yx = -temp/det;
    }
    
    public void E(etomica.space.Tensor t) {E((Tensor2D)t);}
    public void E(etomica.space.Vector u1, etomica.space.Vector u2) {E((Vector2D)u1, (Vector2D)u2);}
    public void PE(etomica.space.Tensor t) {PE((Tensor2D)t);}
    public void PE(etomica.space.Vector u1, etomica.space.Vector u2) {PE((Vector2D)u1, (Vector2D)u2);}
    public void TE(double a) {xx*=a; xy*=a; yx*=a; yy*=a;}
    public void TE(Tensor t){ Tensor2D u = (Tensor2D)t;
	    double txx=xx;double txy=xy;
	    double tyx=yx;double tyy=yy;
	    xx= txx*u.xx+txy*u.yx;   
	    xy= txx*u.xy+txy*u.yy; 
	    yx= tyx*u.xx+tyy*u.yx;   
	    yy= tyx*u.xy+tyy*u.yy; 
    }
    public void E(double[] d) {
        if(d.length != 4) throw new IllegalArgumentException("Array size incorrector for tensor");
        xx = d[0]; xy = d[1]; 
        yx = d[2]; yy = d[3];
    }
    public void assignTo(double[] d) {
        if(d.length != 1) throw new IllegalArgumentException("Array size incorrector for tensor");
        d[0] = xx; d[1] = xy; 
        d[2] = yx; d[3] = yy;
    }

}
