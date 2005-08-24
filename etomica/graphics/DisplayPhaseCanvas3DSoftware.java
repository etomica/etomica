package etomica.graphics;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.TextField;
import java.awt.image.MemoryImageSource;
import java.util.Hashtable;
import java.util.Iterator;

import etomica.Species;
import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.space.Boundary;
import etomica.space.Vector;

    /* History of changes
     * 07/16/02 (DAK) Modified for AtomType.Sphere diameter and radius method to take atom as argument.
     * 09/07/02 (DAK) Added atom filter.
     */

//Class used to define canvas onto which configuration is drawn
public class DisplayPhaseCanvas3DSoftware extends DisplayCanvas {
    private TextField scaleText = new TextField();
    private Font font = new Font("sansserif", Font.PLAIN, 10);
    //  private int annotationHeight = font.getFontMetrics().getHeight();
    private int annotationHeight = 12;
    private static Hashtable ballImages = new Hashtable();
    private final static Color wellColor = new Color(185,185,185, 110);
	/* Atom depth factor. */
	private double atomDepthFactor = 0.33;
	/* Atom sphere factor.*/
	private double atomSphereFactor = 1.;
    private double atomZOffset = .5;

    public Matrix3D amat = new Matrix3D(); // Matrix to do mouse angular rotations.
    public Matrix3D tmat = new Matrix3D(); // Matrix to do translations.
    public Matrix3D zmat = new Matrix3D(); // Matrix to do zooming.
    public Matrix3D mat = new Matrix3D();  // Final matrix for assembly on screen.
    private float prevx, prevy;
    
    private double xfac, xcenter, ycenter, zcenter;

    private int nvert;          // The number of atoms
    private Atom atoms[];       // The entire group of atoms
    private double vert[];      // The verticies of said atoms
    private int tvert[];        // The atom positions transformed to screen space
    private int ZsortMap[];
    private int[] shiftOrigin = new int[2];     //work vector for drawing overflow images

    public void setPrevX(float x) {prevx = x;}
    public void setPrevY(float y) {prevy = y;}
    public float getPrevX() {return(prevx);}
    public float getPrevY() {return(prevy);}
    
    private ColorScheme colorScheme;
    private AtomFilter atomFilter;
        
    public DisplayPhaseCanvas3DSoftware(DisplayPhase _phase) {
        tvert = new int[4];
        scaleText.setVisible(true);
        scaleText.setEditable(false);
        scaleText.setBounds(0,0,100,50);
        displayPhase = _phase;
        colorScheme = displayPhase.getColorScheme();
        atomFilter = displayPhase.getAtomFilter();
    }
    
    public void setAtomFilter(AtomFilter filter) {atomFilter = filter;}
    
    public void initialize() {
        nvert = displayPhase.getPhase().getSpeciesMaster().moleculeCount();
        vert = new double[nvert*3];
        tvert = new int[nvert*3];
        atoms = new Atom[nvert];
        ZsortMap = new int[nvert];
        for (int i = nvert; --i >= 0;)
            ZsortMap[i] = i * 3;
        double xmin = 0., xmax = 0.;
        double ymin = 0., ymax = 0.;
        double zmin = 0., zmax = 0.;
        //temporary commenting
        //this must be put back in for component to work
        /*
        for(Atom a = displayPhase.getPhase().firstAtom(); a!=null; a=a.nextAtom(), k+=3) {
            atoms[k/3] = a;
            vert[k] = a.coord.position().component(0);
            vert[k+1] = a.coord.position().component(1);
            vert[k+2] = a.coord.position().component(2);
            if (vert[k] < xmin) xmin = vert[k];
            if (vert[k] > xmax) xmax = vert[k];
            if (vert[k+1] < ymin) ymin = vert[k+1];
            if (vert[k+1] > ymax) ymax = vert[k+1];
            if (vert[k+2] < zmin) zmin = vert[k+2];
            if (vert[k+2] > zmax) zmax = vert[k+2];
        }
        */
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
    
    private static float[] normalize(float v[]) {
        float len = (float) Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        float v2[] = new float[3];
        if (len == 0.0F) {
            v2[0] = 0.0F;
            v2[1] = 0.0F;
            v2[2] = 0.0F;
        } else {
            v2[0] = v[0] / len;
            v2[1] = v[1] / len;
            v2[2] = v[2] / len;
        }
        return v2;
    }

    private Image SphereSetup(float[] lightSourceVector, java.awt.Color ballColor) {
        float v1[] = new float[3];
        float v2[] = new float[3];
        byte b = 40;
        int i = 2*b + 1;
        int j = -1;
        // Create our own version of an IndexColorModel:
        int model[] = new int[i*i];
        // Normalize the lightsource vector:
        float[] lightsource = normalize(lightSourceVector);
        /*preserve transparency of the color */
        int R2, G2, B2, A2 = ballColor.getAlpha();
        A2 = A2 << 24;
        for (int k1 = -b; k1 <= b; k1++) {
            for (int k2 = -b; k2 <= b; k2++) {
                j++;
                v1[0] = k2;
                v1[1] = k1;
                float len1 = (float) Math.sqrt(k2*k2 + k1*k1);
                if (len1 <=b) {
                    R2 = 0;
                    G2 = 0;
                    B2 = 0;
                    v1[2] = (float)b * (float)Math.cos(Math.asin(len1/b));
                    v1 = normalize(v1);
                    float len2 = (float)Math.abs((double)(v1[0]*lightsource[0] 
                                                          +v1[1]*lightsource[1]
                                                          +v1[2]*lightsource[2]
                                                          ));
                    if (len2 < 0.995f) {
                        R2 = (int)((float)ballColor.getRed()*len2);
                        G2 = (int)((float)ballColor.getGreen()*len2);
                        B2 = (int)((float)ballColor.getBlue()*len2);
                    } else {
                        v2[0] = lightsource[0]+0.0f;
                        v2[1] = lightsource[1]+0.0f;
                        v2[2] = lightsource[2]+1.0f;
                        v2 = normalize(v2);
                        float len3 = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
                        float len4 = 8.0f * len3*len3 - 7.0f;
                        float len5 = 100.0f * len4;
                        len5 = Math.max(len5, 0.0f);                        
                        R2 = (int)((float)(ballColor.getRed() * 155) * len2 + 100.0 + len5);
                        G2 = (int)((float)(ballColor.getGreen() * 155) * len2 + 100.0 + len5);
                        B2 = (int)((float)(ballColor.getBlue() * 155) * len2 + 100.0 + len5);
                        R2 = Math.min(R2, 255);
                        G2 = Math.min(G2, 255);
                        B2 = Math.min(B2, 255);
                    }
                    // Bitwise masking to make model:
                    //model[j] = -16777216 | R2 << 16 | G2 << 8 | B2;
                    model[j] = A2 | R2 << 16 | G2 << 8 | B2;
                } else 
                    model[j] = 0;
            }            
        }
        return createImage(new MemoryImageSource(i, i, model, 0, i));
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
    //        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
            displayPhase.setScale(displayPhase.getScale()*factor);
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
       
	/*
	* Returns the on-screen radius of an atom with the given radius.
	*
	* @param z z position in screen space
	*/
	private double getCircleRadius(double z, double radius) {
		double depth = (z - atomZOffset) / (2.0 * atomZOffset);
	    return ((radius * atomSphereFactor) + (atomDepthFactor * depth) + 7.);
		//double raw = radius * atomSphereFactor;
		//float tmp = atomScreenScale * ((float) raw + atomDepthFactor * depth);
		//double tmp = /*displayPhase.getScale() * displayPhase.getToPixels() */ (raw + atomDepthFactor * depth)+8.;
		//return (tmp > 1.0) ? tmp: 1.0;
	}
	    
    /* Transform all the points in this model */
    void transform() {
        if (nvert <= 0)
            return;
        if (tvert == null || tvert.length < nvert * 3)
            tvert = new int[nvert * 3];
        mat.transform(vert, tvert, nvert);
    }
        
    private void drawAtom(Graphics g, int origin[], Atom a, int sorted) {
        
        if(!atomFilter.accept(a)) return;
        
	    /*  Place the light source for shaded atoms to the upper right of
	    * the atoms and out of the plane. */
        float[] lightSource = {1.0f, -1.0f, 2.0f};
        Vector r = a.coord.position();
        boolean drawWell = false, drawOrientation = false;
        int sigmaP, xP, yP, baseXP, baseYP;
        Image shadedImage;

        if(a.type instanceof AtomTypeOrientedSphere) {
            drawOrientation = true;
        } else if(a.type instanceof AtomTypeWell) {
            drawWell = true;
        }

        Color atomColor = colorScheme.atomColor(a);
        g.setColor(atomColor);
            
        if(a.type instanceof AtomTypeSphere) {
            baseXP = origin[0] + (int)(displayPhase.getToPixels()*(tvert[sorted]+(displayPhase.getPhase().boundary().dimensions().x(0)/2)));
            baseYP = origin[1] + (int)(displayPhase.getToPixels()*(tvert[sorted+1]+(displayPhase.getPhase().boundary().dimensions().x(1)/2)));
            /* Draw the core of the atom */
            sigmaP = (int)(2.0f*getCircleRadius(tvert[sorted+2], ((AtomTypeSphere)a.type).radius(a)));
            //sigmaP = (int)(displayPhase.getToPixels()*2.0f*((AtomType.Sphere)a.type).radius());
            xP = baseXP - (sigmaP>>1);
            yP = baseYP - (sigmaP>>1);
            if (ballImages.containsKey(atomColor)) {
                shadedImage = (Image)ballImages.get(atomColor);
            } else {
                shadedImage = SphereSetup(lightSource, atomColor);
                ballImages.put(atomColor, shadedImage);
            }
            g.drawImage(shadedImage, xP, yP, sigmaP, sigmaP, this);
            /* Draw the surrounding well, if any */
            if(drawWell) {
                g.setColor(wellColor);
                sigmaP = (int)(2.0f*getCircleRadius(tvert[sorted+2], ((AtomTypeWell)a.type).wellRadius()));
                xP = baseXP - (sigmaP>>1);
                yP = baseYP - (sigmaP>>1);
                if (ballImages.containsKey(wellColor)) {
                    shadedImage = (Image)ballImages.get(wellColor);
                } else {
                    shadedImage = SphereSetup(lightSource, wellColor);
                    ballImages.put(wellColor, shadedImage);
                }
                g.drawImage(shadedImage, xP, yP, sigmaP, sigmaP, this);
            }
            /* Draw the orientation line, if any */
            // Needs to be fixed, also needs a true 3d engine to really work.
            /*if(drawOrientation) {
                double theta = ((Space.Coordinate.Angular)a.coordinate()).orientation().angle()[0];
                int dxy = (int)(displayPhase.getToPixels()*getCircleRadius(tvert[sorted+2], ((AtomType.OrientedSphere)a.type).radius()));
                int dx = (int)(dxy*Math.cos(theta));
                int dy = (int)(dxy*Math.sin(theta));
                g.setColor(Color.red);
                xP += dxy; yP += dxy;
                g.drawLine(xP-dx, yP-dy, xP+dx, yP+dy);
            }*/
//            a.type.electroType().draw(g, origin, displayPhase.getToPixels(), r);
        } else { // Not a sphere, wall, or one of their derivatives...
            // Do nothing (how do you draw an object of unkown shape?)
        }
    }
            
    protected boolean computeShiftOrigin(Atom a, Boundary b) {
        if(a.type instanceof AtomTypeSphere) {
            float[][] shifts = b.getOverflowShifts(a.coord.position(),((AtomTypeSphere)a.type).radius(a));  //should instead of radius have a size for all AtomC types
            for(int i=0; i<shifts.length; i++) {
                shiftOrigin[0] = displayPhase.getOrigin()[0] + (int)(displayPhase.getToPixels()*shifts[i][0]);
                shiftOrigin[1] = displayPhase.getOrigin()[1] + (int)(displayPhase.getToPixels()*shifts[i][1]);
            }
            return(true);
        }
        return(false);
    }
      
    /**
    * doPaint is the method that handles the drawing of the phase to the screen.
    * Several variables and conditions affect how the image is drawn.  First,
    * the Unit.Length.Sim class variable <code>TO_PIXELS</code> performs the conversion between simulation
    * length units (Angstroms) and pixels.  The default value is 10 pixels/Angstrom
    * reflecting the default size of the phase (300 pixels by 300 pixels) and the
    * default phase size (30 by 30 A).  
    * The field <code>scale</code> is a multiplicative factor that directly
    * scales up or down the size of the image; this value is adjusted automatically
    * whenever shells of periodic images are drawn, to permit the central image and all 
    * of the specified periodic images to fit in the drawing of the phase.  
    *
    * @param g The graphic object to which the image of the phase is drawn
    * @see Species
    */
    public void doPaint(Graphics g) {
        if(!isVisible() || displayPhase.getPhase() == null) {return;}
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
        displayPhase.computeImageParameters2(w, h);

        //Draw other features if indicated
//        if(drawBoundary>DRAW_BOUNDARY_NONE) {displayPhase.getPhase().boundary().draw(g, displayPhase.getOrigin(), displayPhase.getScale());}

        // Update Atom Positions
        for(int i = 0, k = 0; i < atoms.length; i++, k+=3) {
            vert[k] = atoms[i].coord.position().x(0);
            vert[k+1] = atoms[i].coord.position().x(1);
            vert[k+2] = atoms[i].coord.position().x(2);
            //System.out.println("Atom "+i+"("+vert[k]+", "+vert[k+1]+", "+vert[k+2]+")");
        }

        transform();
            
        /*
        * I use a bubble sort since from one iteration to the next, the sort
        * order is pretty stable, so I just use what I had last time as a
        * "guess" of the sorted order.  With luck, this reduces O(N log N)
        * to O(N)
        */
        boolean flipped;
        int a, b;
        for (int i = nvert - 1; --i >= 0;) {
            flipped = false;
            for (int j = 0; j <= i; j++) {
                a = ZsortMap[j];
                b = ZsortMap[j + 1];
                if (tvert[a + 2] > tvert[b + 2]) {
                    ZsortMap[j + 1] = a;
                    ZsortMap[j] = b;
                    flipped = true;
                }
            }
            if (!flipped)
                break;
        }

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayPhase.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, displayPhase.getOrigin(), displayPhase.getScale());
        }
            
        //Color all atoms according to colorScheme in DisplayPhase
        //if(!initialized)
//        displayPhase.getColorScheme().colorAllAtoms();
            
        //Draw all atoms
        if (nvert <= 0) return;
        for (int i = 0, j; i < nvert; i++) {
            j = ZsortMap[i];
            drawAtom(g, displayPhase.getOrigin(), atoms[j/3], j);
        }
            
        //Draw overflow images if so indicated
        //This needs some work to make more general
        /*if(displayPhase.getDrawOverflow() && displayPhase.parentSimulation().space().D() != 3) {
        /    for(Atom a=displayPhase.getPhase().firstAtom(); a!=null; a=a.nextAtom()) {
                if(computeShiftOrigin(a, boundary)) {
                    drawAtom(g, shiftOrigin, a, 0);
                }
            }
        }*/

        //Draw periodic images if indicated
        /*if(displayPhase.getImageShells() > 0) {
            double[][] origins = displayPhase.getPhase().boundary().imageOrigins(displayPhase.getImageShells());  //more efficient to save rather than recompute each time
            for(int i=0; i<origins.length; i++) {
                g.copyArea(displayPhase.getOrigin()[0],displayPhase.getOrigin()[1],displayPhase.getDrawSize()[0],displayPhase.getDrawSize()[1],(int)(displayPhase.getToPixels()*origins[i][0]),(int)(displayPhase.getToPixels()*origins[i][1]));
            }
        }*/
        //Draw bar showing scale if indicated
        if(writeScale) {
            g.setColor(Color.lightGray);
            g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
            g.setColor(Color.black);
            g.setFont(font);
            g.drawString("Scale: "+Integer.toString((int)(100*displayPhase.getScale()))+"%", 0, getSize().height-3);
        }
        //System.gc();
    }

        /** A fairly conventional 3D matrix object that can transform sets of
        3D points and perform a variety of manipulations on the transform */
    public static class Matrix3D implements java.io.Serializable {
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
    
}  //end of DisplayPhase.Canvas
