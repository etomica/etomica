/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;


import etomica.math.function.IFunction;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.Debug;

public class Tensor3D implements Tensor, java.io.Serializable {

    // xy is row x column y
    protected double xx, xy, xz, yx, yy, yz, zx, zy, zz;
    private static final long serialVersionUID = 1L;

    /**
     * Default constructor sets all elements to zero.
     */
    public Tensor3D () {
        xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;
    }
    
    
    /**
     * Constructs tensor with elements set by the given array.  Elements
     * are interpreted in order as xx, xy, xz, yx, yy, yz, zx, zy, zz.
     */
    public Tensor3D (double[][] d) {
        this.E(d);
    }
    
    /**
     * Returns tensor in a new array, writing elements row-wise, i.e.: xx, xy, xz, yx, yy, yz, zx, zy, zz;
     */
    public double[] toArray() {
        return new double[] {xx, xy, xz, yx, yy, yz, zx, zy, zz};
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
        return ( i==0 ) ? ( (j==0) ? xx : ( j==1 ? xy : xz ) ) : ( (i==1) ? ( (j==0) ? yx : ( (j==1) ? yy : yz ) ) : ( (j==0) ? zx : ((j==1) ? zy : zz)));
    }
    
    public int D() {return 3;}
    
    public void setComponent(int i, int j, double d) {
        if (i==0) {if (j==0) {xx = d;} else if (j==1) {xy = d;} else xz = d;}
        else if (i==1) {if (j==0) {yx = d;} else if (j==1) {yy=d;} else yz = d;}
        else {if (j==0) {zx = d;} else if (j==1) {zy = d;} else zz = d;}
    }

    public void E(Tensor u) {
        Tensor3D t = (Tensor3D)u;
        xx=t.xx; xy=t.xy; xz=t.xz;
        yx=t.yx; yy=t.yy; yz=t.yz;
        zx=t.zx; zy=t.zy; zz=t.zz;
    }
    
    public void E(Vector[] v) {
        if(v.length != 3) {
            throw new IllegalArgumentException("Tensor requires 3 vectors to set its values");
        }
        xx = ((Vector3D)v[0]).x; xy = ((Vector3D)v[1]).x; xz = ((Vector3D)v[2]).x;
        yx = ((Vector3D)v[0]).y; yy = ((Vector3D)v[1]).y; yz = ((Vector3D)v[2]).y;
        zx = ((Vector3D)v[0]).z; zy = ((Vector3D)v[1]).z; zz = ((Vector3D)v[2]).z;
    }

    public void E(double[] d) {
        if(Debug.ON && d.length != 9) throw new IllegalArgumentException("Array size incorrect for tensor");
        xx = d[0]; xy = d[1]; xz = d[2];
        yx = d[3]; yy = d[4]; yz = d[5];
        zx = d[6]; zy = d[7]; zz = d[8];
    }

    public void E(double[][] d) {
        if(Debug.ON && d.length != 3) throw new IllegalArgumentException("Array size incorrect for tensor");
        xx = d[0][0]; xy = d[0][1]; xz = d[0][2];
        yx = d[1][0]; yy = d[1][1]; yz = d[1][2];
        zx = d[2][0]; zy = d[2][1]; zz = d[2][2];
    }


    public void diagE(Vector v) {
        this.E(0.0);
        Vector3D v3 = (Vector3D)v;
        xx = v3.x;
        yy = v3.y;
        zz = v3.z;
    }
    
    public void assignTo(Vector[] v) {
        if (v.length != 3) {
            throw new IllegalArgumentException("Tensor requires 3 vector for assignment");
        }
        double[][] a = new double[3][3];
        a[0][0] = xx;
        a[1][0] = xy;
        a[2][0] = xz;
        a[0][1] = yx;
        a[1][1] = yy;
        a[2][1] = yz;
        a[0][2] = zx;
        a[1][2] = zy;
        a[2][2] = zz;
        for (int i = 0; i < 3; i++) {
            v[i].E(a[i]);
        }
    }

    public void Ev1v2(Vector v1, Vector v2) {
        xx=v1.x()*v2.x(); xy=v1.x()*v2.y(); xz=v1.x()*v2.z();
        yx=v1.y()*v2.x(); yy=v1.y()*v2.y(); yz=v1.y()*v2.z();
        zx=v1.z()*v2.x(); zy=v1.z()*v2.y(); zz=v1.z()*v2.z();
    }
    
    public void E(double a) {
        xx=xy=xz=yx=yy=yz=zx=zy=zz=a;
    }
    
    public void PE(double a) {
        xx+=a; xy+=a; xz+=a;
        yx+=a; yy+=a; yz+=a;
        zx+=a; zy+=a; zz+=a;
    }
    
    public void PE(Tensor u) {
        Tensor3D t = (Tensor3D)u;
        xx+=t.xx; xy+=t.xy; xz+=t.xz;
        yx+=t.yx; yy+=t.yy; yz+=t.yz;
        zx+=t.zx; zy+=t.zy; zz+=t.zz;
    }
    
    public void PE(int i, int j, double d) {
        if (i==0) {if (j==0) {xx += d;} else if (j==1) {xy += d;} else xz += d;}
        else if (i==1) {if (j==0) {yx += d;} else if (j==1) {yy += d;} else yz += d;}
        else {if (j==0) {zx += d;} else if (j==1) {zy += d;} else zz += d;}
    }
    
    public void ME(Tensor u) {
        Tensor3D t = (Tensor3D)u;
        xx-=t.xx; xy-=t.xy; xz-=t.xz;
        yx-=t.yx; yy-=t.yy; yz-=t.yz;
        zx-=t.zx; zy-=t.zy; zz-=t.zz;
    }
    
    public double trace() {
        return xx+yy+zz;
    }
    
    public void transpose() { 
        	double temp = 0.0;
        	temp = xy; xy = yx; yx = temp;
        	temp = xz; xz = zx; zx = temp;
        	temp = zy; zy = yz; yz = temp;    	
    }
    
    public double determinant() {
        return xx*yy*zz-xx*yz*zy-xy*yx*zz+xz*yx*zy+xy*yz*zx-xz*yy*zx;
    }

    public void invert() {
        double txx=xx;double txy=xy;double txz=xz;
        double tyx=yx;double tyy=yy;double tyz=yz;
        double tzx=zx;double tzy=zy;double tzz=zz;
	    double det = determinant();
        xx= (tyy*tzz-tyz*tzy)/det; 
        xy= -(txy*tzz-txz*tzy)/det;
        xz= (txy*tyz-txz*tyy)/det;
        yx= -(tyx*tzz-tyz*tzx)/det; 
        yy= (txx*tzz-txz*tzx)/det;
        yz= -(txx*tyz-txz*tyx)/det;
        zx= (tyx*tzy-tyy*tzx)/det;
        zy= -(txx*tzy-txy*tzx)/det;
        zz= (txx*tyy-txy*tyx)/det;                              
    }
    
    public void PEv1v2(Vector v1, Vector v2) {
        Vector3D u1 = (Vector3D)v1;
        Vector3D u2 = (Vector3D)v2;
        xx+=u1.x*u2.x; xy+=u1.x*u2.y; xz+=u1.x*u2.z;
        yx+=u1.y*u2.x; yy+=u1.y*u2.y; yz+=u1.y*u2.z;
        zx+=u1.z*u2.x; zy+=u1.z*u2.y; zz+=u1.z*u2.z;
    }
    
    public void MEv1v2(Vector v1, Vector v2) {
        Vector3D u1 = (Vector3D)v1;
        Vector3D u2 = (Vector3D)v2;
        xx-=u1.x*u2.x; xy-=u1.x*u2.y; xz-=u1.x*u2.z;
        yx-=u1.y*u2.x; yy-=u1.y*u2.y; yz-=u1.y*u2.z;
        zx-=u1.z*u2.x; zy-=u1.z*u2.y; zz-=u1.z*u2.z;
    }
    
    public void PEa1Tt1(double a1, Tensor t1) {
        xx += a1*((Tensor3D)t1).xx;
        xy += a1*((Tensor3D)t1).xy;
        xz += a1*((Tensor3D)t1).xz;
        yx += a1*((Tensor3D)t1).yx;
        yy += a1*((Tensor3D)t1).yy;
        yz += a1*((Tensor3D)t1).yz;
        zx += a1*((Tensor3D)t1).zx;
        zy += a1*((Tensor3D)t1).zy;
        zz += a1*((Tensor3D)t1).zz;
    }

    public void TE(double a) {
        xx*=a; xy*=a; xz*=a; 
        yx*=a; yy*=a; yz*=a; 
        zx*=a; zy*=a; zz*=a;
    }
    
    public void TE(Tensor t) { Tensor3D u = (Tensor3D)t;
        double txx=xx;double txy=xy;double txz=xz;
        double tyx=yx;double tyy=yy;double tyz=yz;
        double tzx=zx;double tzy=zy;double tzz=zz;
        xx= txx*u.xx+txy*u.yx+txz*u.zx;   
        xy= txx*u.xy+txy*u.yy+txz*u.zy; 
        xz= txx*u.xz+txy*u.yz+txz*u.zz;
        yx= tyx*u.xx+tyy*u.yx+tyz*u.zx;   
        yy= tyx*u.xy+tyy*u.yy+tyz*u.zy; 
        yz= tyx*u.xz+tyy*u.yz+tyz*u.zz;
        zx= tzx*u.xx+tzy*u.yx+tzz*u.zx;   
        zy= tzx*u.xy+tzy*u.yy+tzz*u.zy; 
        zz= tzx*u.xz+tzy*u.yz+tzz*u.zz;                                         
    }
    
    public void DE(Tensor t) {
        Tensor3D u = (Tensor3D)t;
        xx /= u.xx; xy /= u.xy; xz /= u.xz;
        yx /= u.yx; yy /= u.yy; yz /= u.yz;
        zx /= u.zx; zy /= u.zy; zz /= u.zz;
    }
    
    public void assignTo(double[] d) {
        if(Debug.ON && d.length != 9) throw new IllegalArgumentException("Array size incorrect for tensor");
        d[0] = xx; d[1] = xy; d[2] = xz; 
        d[3] = yx; d[4] = yy; d[5] = yz;
        d[6] = zx; d[7] = zy; d[8] = zz;
    }
    
    public void assignTo(double[][] d) {
        d[0][0] = xx; d[0][1] = xy; d[0][2] = xz;
        d[1][0] = yx; d[1][1] = yy; d[1][2] = yz;
        d[2][0] = zx; d[2][1] = zy; d[2][2] = zz;
    }

    public boolean isNaN() {
        return Double.isNaN(xx) || Double.isNaN(xy) || Double.isNaN(xz)
            || Double.isNaN(yx) || Double.isNaN(yy) || Double.isNaN(yz)
            || Double.isNaN(zx) || Double.isNaN(zy) || Double.isNaN(zz);
    }
    
    public void map(IFunction f) {
        xx = f.f(xx); xy = f.f(xy); xz = f.f(xz);
        yx = f.f(yx); yy = f.f(yy); yz = f.f(yz);
        zx = f.f(zx); zy = f.f(zy); zz = f.f(zz);
    }

    public void transform(Vector v) {
        Vector3D v3D = (Vector3D) v;
        double x1 = xx * v3D.x + xy * v3D.y + xz * v3D.z;
        double y1 = yx * v3D.x + yy * v3D.y + yz * v3D.z;
        v3D.z = zx * v3D.x + zy * v3D.y + zz * v3D.z;
        v3D.x = x1;
        v3D.y = y1;
    }

    public boolean equals(Tensor t) {
        Tensor3D A = (Tensor3D)t;
        return A.xx==xx && A.xy==xy && A.xz==xz && A.yx==yx && A.yy==yy && A.yz==yz && A.zx==zx && A.zy==zy && A.zz==zz;
    }

    public String toString() {
        return "("+xx+", "+xy+", "+xz+")\n("+yx+", "+yy+", "+yz+")"+"\n("+zx+", "+zy+", "+zz+")";
    }

}
