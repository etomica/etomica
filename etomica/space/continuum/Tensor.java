package etomica.space.continuum;
import etomica.Space;

public class Tensor implements Space.Tensor {
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
    public static final Tensor ORIGIN = new Tensor();
    public static final Tensor WORK = new Tensor();
    public Tensor () {xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;}
    public Tensor (double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
        this.xx=xx; this.xy=xy; this.xz=xz; this.yx=yx; this.yy=yy; this.yz=yz; this.zx=zx; this.zy=zy; this.zz=zz;
    }
    public double component(int i, int j) {
        return ( i==0 ) ? ( (j==0) ? xx : ( j==1 ? xy : xz ) ) : ( (i==1) ? ( (j==0) ? yx : ( (j==1) ? yy : yz ) ) : ( (j==0) ? zx : ((j==1) ? zy : zz)));
    }
    public int length() {return 3;}
    public void setComponent(int i, int j, double d) {
        if (i==0) {if (j==0) {xx = d;} else if (j==1) {xy = d;} else xz = d;}
        else if (i==1) {if (j==0) {yx = d;} else if (j==1) {yy=d;} else yz = d;}
        else {if (j==0) {zx = d;} else if (j==1) {zy = d;} else zz = d;}
    }
    public void E(Tensor t) {xx=t.xx; xy=t.xy; xz=t.xz; yx=t.yx; yy=t.yy; yz=t.yz; zx=t.zx; zy=t.zy; zz=t.zz;}
    public void E(Vector u1, Vector u2) {
        double[] x1 = u1.toArray();
        double[] x2 = u2.toArray();
        xx=x1[0]*x2[0]; xy=x1[0]*x2[1]; xz=x1[0]*x2[2]; 
        yx=x1[1]*x2[0]; yy=x1[1]*x2[1]; yz=x1[1]*x2[2]; 
        zx=x1[2]*x2[0]; zy=x1[2]*x2[1]; zz=x1[2]*x2[2];
    }
    public void E(double a) {xx=xy=xz=yx=yy=yz=zx=zy=zz=a;}
    public void PE(Tensor t) {xx+=t.xx; xy+=t.xy; xz+=t.xz; yx+=t.yx; yy+=t.yy; yz+=t.yz; zx+=t.zx; zy+=t.zy; zz+=t.zz;}
    public void PE(int i, int j, double d) {
        if (i==0) {if (j==0) {xx += d;} else if (j==1) {xy += d;} else xz += d;}
        else if (i==1) {if (j==0) {yx += d;} else if (j==1) {yy += d;} else yz += d;}
        else {if (j==0) {zx += d;} else if (j==1) {zy += d;} else zz += d;}
    }
    public double trace() {return xx+yy+zz;}
    public void E(Space.Tensor t) {E((Tensor)t);}
    public void E(Space.Vector u1, Space.Vector u2) {E((Vector)u1, (Vector)u2);}
    public void PE(Space.Tensor t) {PE((Tensor) t);}
    public void PE(Space.Vector u1, Space.Vector u2) {PE((Vector)u1, (Vector)u2);}
    public void TE(double a) {xx*=a; xy*=a; xz*=a; yx*=a; yy*=a; yz*=a; zx*=a; zy*=a; zz*=a;}
}
    
