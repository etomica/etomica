package etomica.graphics;
import etomica.*;

import java.awt.*;
import java.awt.image.MemoryImageSource;

import java.util.Hashtable;
//java2 replacement
import etomica.utility.java2.Iterator;

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
        nvert = 0;
        for(Iterator iter=displayPhase.simulation().speciesList().iterator(); iter.hasNext(); ) {
            Species obj = (Species)iter.next();
            nvert += obj.getNMolecules();
        }
        vert = new double[nvert*3];
        tvert = new int[nvert*3];
        atoms = new Atom[nvert];
        ZsortMap = new int[nvert];
        for (int i = nvert; --i >= 0;)
            ZsortMap[i] = i * 3;
        int k = 0;
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
        Space.Vector r = a.coord.position();
        boolean drawWell = false, drawOrientation = false;
        int sigmaP, xP, yP, baseXP, baseYP;
        Image shadedImage;

        if(a.type instanceof AtomType.OrientedSphere) {
            drawOrientation = true;
        } else if(a.type instanceof AtomType.Well) {
            drawWell = true;
        }

        Color atomColor = colorScheme.atomColor(a);
        g.setColor(atomColor);
            
        if(a.type instanceof AtomType.Sphere) {
            baseXP = origin[0] + (int)(displayPhase.getToPixels()*(tvert[sorted]+(displayPhase.getPhase().boundary().dimensions().x(0)/2)));
            baseYP = origin[1] + (int)(displayPhase.getToPixels()*(tvert[sorted+1]+(displayPhase.getPhase().boundary().dimensions().x(1)/2)));
            /* Draw the core of the atom */
            sigmaP = (int)(2.0f*getCircleRadius(tvert[sorted+2], ((AtomType.Sphere)a.type).radius(a)));
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
                sigmaP = (int)(2.0f*getCircleRadius(tvert[sorted+2], ((AtomType.Well)a.type).wellRadius()));
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
        } else if(a.type instanceof AtomType.Wall) {
            xP = origin[0] + (int)(displayPhase.getToPixels()*r.x(0));
            yP = origin[1] + (int)(displayPhase.getToPixels()*r.x(1));
            int t = Math.max(1,(int)((double)((AtomType.Wall)a.type).getThickness()*(double)displayPhase.getToPixels()/(double)etomica.units.BaseUnit.Length.Sim.TO_PIXELS));
            if(!(((AtomType.Wall)a.type).isHorizontal() || ((AtomType.Wall)a.type).isVertical())) {  //not horizontal or vertical; draw line
                int x1 = xP + (int)(displayPhase.getToPixels()*((AtomType.Wall)a.type).getLength()*((AtomType.Wall)a.type).getCosX());
                int y1 = yP + (int)(displayPhase.getToPixels()*((AtomType.Wall)a.type).getLength()*((AtomType.Wall)a.type).getSinX());
                g.drawLine(xP, yP, x1, y1);
            }
            else if(((AtomType.Wall)a.type).isLongWall()) {
                java.awt.Rectangle rect = g.getClipBounds();
                //int wP = vertical ? t : (int)(toPixels*atom.parentPhase().boundary().dimensions().x(1));
                //int hP = horizontal ? t : (int)(toPixels*atom.parentPhase().boundary().dimensions().x(0));
                int wP = ((AtomType.Wall)a.type).isVertical() ? t : Integer.MAX_VALUE;
                int hP = ((AtomType.Wall)a.type).isHorizontal() ? t : Integer.MAX_VALUE;
                //int X = vertical ? xP : origin[0];
                //int Y = horizontal ? yP : origin[1];
                int X = ((AtomType.Wall)a.type).isVertical() ? xP : 0;
                int Y = ((AtomType.Wall)a.type).isHorizontal() ? yP : 0;
                g.fillRect(X,Y,wP,hP);
            }   
            else {                           //horizontal or vertical; draw box
                int wP = ((AtomType.Wall)a.type).isVertical() ? t : (int)(displayPhase.getToPixels()*((AtomType.Wall)a.type).getLength());
                int hP = ((AtomType.Wall)a.type).isHorizontal() ? t : (int)(displayPhase.getToPixels()*((AtomType.Wall)a.type).getLength());
                g.fillRect(xP,yP,wP,hP);
            }
        } else { // Not a sphere, wall, or one of their derivatives...
            // Do nothing (how do you draw an object of unkown shape?)
        }
    }
            
    protected boolean computeShiftOrigin(Atom a, Space.Boundary b) {
        if(a.type instanceof AtomType.Sphere) {
            float[][] shifts = b.getOverflowShifts(a.coord.position(),((AtomType.Sphere)a.type).radius(a));  //should instead of radius have a size for all AtomC types
            for(int i=0; i<shifts.length; i++) {
                shiftOrigin[0] = displayPhase.getOrigin()[0] + (int)(displayPhase.getToPixels()*shifts[i][0]);
                shiftOrigin[1] = displayPhase.getOrigin()[1] + (int)(displayPhase.getToPixels()*shifts[i][1]);
            }
            return(true);
        } else {
            return(false);
        }
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
        Space.Boundary boundary = displayPhase.getPhase().boundary();
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
}  //end of DisplayPhase.Canvas
