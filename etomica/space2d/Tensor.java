package simulate.space2D;
import simulate.Space;

public final class Tensor extends Space.Tensor {
    public int length() {return 2;}
    double xx, xy, yx, yy;
    public static final Tensor ZERO = new Tensor();
    public static final Tensor WORK = new Tensor();  //anything using WORK is not thread-safe
    public Tensor () {xx = xy = yx = yy = 0.0;}
    public Tensor (double xx, double xy, double yx, double yy) {this.xx = xx; this.xy = xy; this.yx = yx; this.yy = yy;}

    public double component(int i, int j) {
        return (i==0) ? ( (j==0) ? xx : xy ) : ( (j==0) ? yx : yy );}
    public void setComponent(int i, int j, double d) {
        if(i==0) {if(j==0) xx=d; else xy=d;}
        else     {if(j==0) yx=d; else yy=d;}
    }
    public void E(Tensor t) {xx=t.xx; xy=t.xy; yx=t.yx; yy=t.yy;}
    public void E(Vector u1, Vector u2) {xx=u1.x*u2.x; xy=u1.x*u2.y; yx=u1.y*u2.x; yy=u1.y*u2.y;}
    public void E(double a) {xx = xy = yx = yy = a;}
    public void PE(Tensor t) {xx+=t.xx; xy+=t.xy; yx+=t.yx; yy+=t.yy;}
    public void PE(int i, int j, double a) {
        if(i==0) {if(j==0) xx+=a; else xy+=a;}
        else     {if(j==0) yx+=a; else yy+=a;}
    }
    public void PE(Vector u1, Vector u2) {xx+=u1.x*u2.x; xy+=u1.x*u2.y; yx+=u1.y*u2.x; yy+=u1.y*u2.y;}
    public double trace() {return xx + yy;}
        
    public void E(Space.Tensor t) {E((Tensor)t);}
    public void E(Space.Vector u1, Space.Vector u2) {E((Vector)u1, (Vector)u2);}
    public void PE(Space.Tensor t) {PE((Tensor)t);}
    public void PE(Space.Vector u1, Space.Vector u2) {PE((Vector)u1, (Vector)u2);}
    public void TE(double a) {xx*=a; xy*=a; yx*=a; yy*=a;}
}

