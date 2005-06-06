package etomica.graphics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import etomica.atom.AtomFilter;
import etomica.exception.MethodNotImplementedException;

//Class used to define canvas onto which configuration is drawn
public class DisplayPolytopeCanvas3DSoftware extends DisplayCanvas {
    public Matrix3D amat = new Matrix3D(); // Matrix to do mouse angular rotations.
    public Matrix3D tmat = new Matrix3D(); // Matrix to do translations.
    public Matrix3D zmat = new Matrix3D(); // Matrix to do zooming.
    public Matrix3D mat = new Matrix3D();  // Final matrix for assembly on screen.
    private float prevx, prevy;
    
    private double xfac, xcenter, ycenter, zcenter;

    public void setPrevX(float x) {prevx = x;}
    public void setPrevY(float y) {prevy = y;}
    public float getPrevX() {return(prevx);}
    public float getPrevY() {return(prevy);}
    
    private DisplayPolytope displayPolytope;
        
    public DisplayPolytopeCanvas3DSoftware(DisplayPolytope _polytope) {
        displayPolytope = _polytope;
    }
    
    public void setAtomFilter(AtomFilter filter) {throw new MethodNotImplementedException("None shall pass");}
    
    public void initialize() {
        double xmin = 0., xmax = 0.;
        double ymin = 0., ymax = 0.;
        double zmin = 0., zmax = 0.;
        double xw = xmax - xmin;
        double yw = ymax - ymin;
        double zw = zmax - zmin;
        if (yw > xw)
            xw = yw;
        if (zw > xw)
            xw = zw;
        double f1 = getSize().width / xw;
        double f2 = getSize().height / xw;
        xfac = 0.7 * (f1 < f2 ? f1 : f2);// * scalefudge;
        xcenter = -(xmin + xmax) / 2;
        ycenter = -(ymin + ymax) / 2;
        zcenter = -(zmin + zmax) / 2;
        mat.unit();
        mat.translate(xcenter, ycenter, zcenter);
        //mat.xrot(-5);
        //mat.yrot(-10);
        //mat.scale(w*xfac, -h*xfac, xfac / getSize().width);
        mat.scale(xfac, -xfac,  16* xfac / getSize().width);
    }
    
    /**
        * Sets the size of the display to a new value and scales the image so that
        * the phase fits in the canvas in the same proportion as before.
        */
    public void scaleSetSize(int width, int height) {
        if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
            double ratio1 = (double)width/(double)getBounds().width;
            double ratio2 = (double)height/(double)getBounds().height;
            double factor = Math.min(ratio1, ratio2);
            displayPolytope.setScale(displayPolytope.getScale()*factor);
            setSize(width, height);
        }
    }
            
          
    //Override superclass methods for changing size so that scale is reset with any size change  
    // this setBounds is ultimately called by all other setSize, setBounds methods
    public void setBounds(int x, int y, int width, int height) {
        if(width == 0 || height == 0) return;
        super.setBounds(x,y,width,height);
        createOffScreen(width,height);
    }

    /**
     * doPaint is the method that handles the drawing of the phase to the screen.
     * Several variables and conditions affect how the image is drawn.  First,
     * the Unit.Length.Sim class variable <code>TO_PIXELS</code> performs the conversion 
     * between polytope dimensions and pixels.  The default value is 10 pixels/unit length
     * reflecting the default size of the phase (300 pixels by 300 pixels) and the
     * default polytope size (30 by 30).
     *
     * @param g The graphic object to which the image of the phase is drawn
     */
    public void doPaint(Graphics g) {
        if(!isVisible() || displayPolytope.getPolytope() == null) {return;}
        int w = getSize().width;
        int h = getSize().height;
            
        String vers = System.getProperty("java.version");
        if (vers.compareTo("1.2") >= 0) {
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                RenderingHints.VALUE_ANTIALIAS_OFF);
            g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                RenderingHints.VALUE_RENDER_SPEED);
            g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, 
                RenderingHints.VALUE_COLOR_RENDER_SPEED);
            if(quality>=DRAW_QUALITY_HIGH) {
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
                    RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setRenderingHint(RenderingHints.KEY_RENDERING, 
                    RenderingHints.VALUE_RENDER_QUALITY);
                g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, 
                    RenderingHints.VALUE_COLOR_RENDER_QUALITY);
            }
        }
            
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        displayPolytope.computeImageParameters2(w, h);

        //TODO actually draw the polytope!!!!
    }

    /**
     * A fairly conventional 3D matrix object that can transform sets of
     * 3D points and perform a variety of manipulations on the transform
     */
    public static class Matrix3D {
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
      
      /** PhaseTranslate the origin */
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
    
        double Nxx = xx * ct + zx * st;
        double Nxy = xy * ct + zy * st;
        double Nxz = xz * ct + zz * st;
        double Nxo = xo * ct + zo * st;
    
        double Nzx = zx * ct - xx * st;
        double Nzy = zy * ct - xy * st;
        double Nzz = zz * ct - xz * st;
        double Nzo = zo * ct - xo * st;
    
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
    
        double Nyx = yx * ct + zx * st;
        double Nyy = yy * ct + zy * st;
        double Nyz = yz * ct + zz * st;
        double Nyo = yo * ct + zo * st;
    
        double Nzx = zx * ct - yx * st;
        double Nzy = zy * ct - yy * st;
        double Nzz = zz * ct - yz * st;
        double Nzo = zo * ct - yo * st;
    
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
    
        double Nyx = yx * ct + xx * st;
        double Nyy = yy * ct + xy * st;
        double Nyz = yz * ct + xz * st;
        double Nyo = yo * ct + xo * st;
    
        double Nxx = xx * ct - yx * st;
        double Nxy = xy * ct - yy * st;
        double Nxz = xz * ct - yz * st;
        double Nxo = xo * ct - yo * st;
    
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
    
}  //end of DisplayPhase.Canvas
