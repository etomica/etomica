/*
 * @(#)Matrix3D.java	1.2 96/12/06
 *
 * Copyright (c) 1994-1996 Sun Microsystems, Inc. All Rights Reserved.
 *
 * Sun grants you ("Licensee") a non-exclusive, royalty free, license to use,
 * modify and redistribute this software in source and binary code form,
 * provided that i) this copyright notice and license appear on all copies of
 * the software; and ii) Licensee does not utilize the software in a manner
 * which is disparaging to Sun.
 *
 * This software is provided "AS IS," without a warranty of any kind. ALL
 * EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY
 * IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR
 * NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL NOT BE
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
 * OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER
 * CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE USE OF
 * OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line control of
 * aircraft, air traffic, aircraft navigation or aircraft communications; or in
 * the design, construction, operation or maintenance of any nuclear
 * facility. Licensee represents and warrants that it will not use or
 * redistribute the Software for such purposes.
 */

package etomica;

/** A fairly conventional 3D matrix object that can transform sets of
    3D points and perform a variety of manipulations on the transform */
public class Matrix3D {
  double xx, xy, xz, xo;
  double yx, yy, yz, yo;
  double zx, zy, zz, zo;
  
  /** Create a new unit matrix */
  public Matrix3D () {
    xx = 1.0f;
    yy = 1.0f;
    zz = 1.0f;
  }
  
  /** Scale by f in all dimensions */
  public void scale(double f) {
    xx *= f;
    xy *= f;
    xz *= f;
    xo *= f;
    yx *= f;
    yy *= f;
    yz *= f;
    yo *= f;
    zx *= f;
    zy *= f;
    zz *= f;
    zo *= f;
  }
  
  /** Scale along each axis independently */
  public void scale(double xf, double yf, double zf) {
    xx *= xf;
    xy *= xf;
    xz *= xf;
    xo *= xf;
    yx *= yf;
    yy *= yf;
    yz *= yf;
    yo *= yf;
    zx *= zf;
    zy *= zf;
    zz *= zf;
    zo *= zf;
  }
  
  /** Translate the origin */
  public void translate(double x, double y, double z) {
    xo += x;
    yo += y;
    zo += z;
  }
  
  /** rotate theta degrees about the y axis */
  public void yrot(double theta) {
    theta *= (Math.PI / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nxx = (double) (xx * ct + zx * st);
    double Nxy = (double) (xy * ct + zy * st);
    double Nxz = (double) (xz * ct + zz * st);
    double Nxo = (double) (xo * ct + zo * st);

    double Nzx = (double) (zx * ct - xx * st);
    double Nzy = (double) (zy * ct - xy * st);
    double Nzz = (double) (zz * ct - xz * st);
    double Nzo = (double) (zo * ct - xo * st);

    xo = Nxo;
    xx = Nxx;
    xy = Nxy;
    xz = Nxz;
    zo = Nzo;
    zx = Nzx;
    zy = Nzy;
    zz = Nzz;
  }
  
  /** rotate theta degrees about the x axis */
  public void xrot(double theta) {
    theta *= (Math.PI / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (double) (yx * ct + zx * st);
    double Nyy = (double) (yy * ct + zy * st);
    double Nyz = (double) (yz * ct + zz * st);
    double Nyo = (double) (yo * ct + zo * st);

    double Nzx = (double) (zx * ct - yx * st);
    double Nzy = (double) (zy * ct - yy * st);
    double Nzz = (double) (zz * ct - yz * st);
    double Nzo = (double) (zo * ct - yo * st);

    yo = Nyo;
    yx = Nyx;
    yy = Nyy;
    yz = Nyz;
    zo = Nzo;
    zx = Nzx;
    zy = Nzy;
    zz = Nzz;
  }
  
  /** rotate theta degrees about the z axis */
  public void zrot(double theta) {
    theta *= (Math.PI / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (double) (yx * ct + xx * st);
    double Nyy = (double) (yy * ct + xy * st);
    double Nyz = (double) (yz * ct + xz * st);
    double Nyo = (double) (yo * ct + xo * st);

    double Nxx = (double) (xx * ct - yx * st);
    double Nxy = (double) (xy * ct - yy * st);
    double Nxz = (double) (xz * ct - yz * st);
    double Nxo = (double) (xo * ct - yo * st);

    yo = Nyo;
    yx = Nyx;
    yy = Nyy;
    yz = Nyz;
    xo = Nxo;
    xx = Nxx;
    xy = Nxy;
    xz = Nxz;
  }
  
  /** Multiply this matrix by a second: M = M*R */
  public void mult(Matrix3D rhs) {
    double lxx = xx * rhs.xx + yx * rhs.xy + zx * rhs.xz;
    double lxy = xy * rhs.xx + yy * rhs.xy + zy * rhs.xz;
    double lxz = xz * rhs.xx + yz * rhs.xy + zz * rhs.xz;
    double lxo = xo * rhs.xx + yo * rhs.xy + zo * rhs.xz + rhs.xo;

    double lyx = xx * rhs.yx + yx * rhs.yy + zx * rhs.yz;
    double lyy = xy * rhs.yx + yy * rhs.yy + zy * rhs.yz;
    double lyz = xz * rhs.yx + yz * rhs.yy + zz * rhs.yz;
    double lyo = xo * rhs.yx + yo * rhs.yy + zo * rhs.yz + rhs.yo;

    double lzx = xx * rhs.zx + yx * rhs.zy + zx * rhs.zz;
    double lzy = xy * rhs.zx + yy * rhs.zy + zy * rhs.zz;
    double lzz = xz * rhs.zx + yz * rhs.zy + zz * rhs.zz;
    double lzo = xo * rhs.zx + yo * rhs.zy + zo * rhs.zz + rhs.zo;

    xx = lxx;
    xy = lxy;
    xz = lxz;
    xo = lxo;

    yx = lyx;
    yy = lyy;
    yz = lyz;
    yo = lyo;

    zx = lzx;
    zy = lzy;
    zz = lzz;
    zo = lzo;
  }

  /** Reinitialize to the unit matrix */
  public void unit() {
    xo = 0;
    xx = 1;
    xy = 0;
    xz = 0;
    yo = 0;
    yx = 0;
    yy = 1;
    yz = 0;
    zo = 0;
    zx = 0;
    zy = 0;
    zz = 1;
  }
  
  /** Transform nvert points from v into tv.  v contains the input
  coordinates in floating point.  Three successive entries in
  the array constitute a point.  tv ends up holding the transformed
  points as integers; three successive entries per point */
  public void transform(double v[], int tv[], int nvert) {
    double lxx = xx, lxy = xy, lxz = xz, lxo = xo;
    double lyx = yx, lyy = yy, lyz = yz, lyo = yo;
    double lzx = zx, lzy = zy, lzz = zz, lzo = zo;
    for (int i = nvert * 3; (i -= 3) >= 0;) {
      double x = v[i];
      double y = v[i + 1];
      double z = v[i + 2];
      tv[i    ] = (int) (x * lxx + y * lxy + z * lxz + lxo);
      tv[i + 1] = (int) (x * lyx + y * lyy + z * lyz + lyo);
      tv[i + 2] = (int) (x * lzx + y * lzy + z * lzz + lzo);
    }
  }
      
  public String toString() {
    return ("[" + xo + "," + xx + "," + xy + "," + xz + ";"
    + yo + "," + yx + "," + yy + "," + yz + ";"
    + zo + "," + zx + "," + zy + "," + zz + "]");
  }
}
