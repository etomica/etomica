//http://www.lanw.com/training/interop/securityurls.htm

//http://www.geocities.com/SiliconValley/Park/5625/opengl/
//http://www.azillionmonkeys.com/windoze/OpenGLvsDirect3D.html

//http://romka.demonews.com/opengl/demos/other_eng.htm
//http://www.oglchallenge.com/

//http://www.sgi.com/software/opengl/advanced97/notes/node196.html#stencilsection
//http://ask.ii.uib.no/ebt-bin/nph-dweb/dynaweb/SGI_Developer/OpenGL_RM
package etomica.graphics;
import java.awt.Color;
import java.util.Iterator;

import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polyhedron;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import gl4java.utils.glut.GLUTEnum;

    /* History of changes
     * 07/16/02 (DAK) Modified for AtomType.Sphere diameter and radius method to take atom as argument.
     * 09/07/02 (DAK) added atomFilter
     * 09/xx/02 (DAK) added display of plane
     * 09/27/02 (DAK) set zoom based on size of phase boundary (in init method)
     * 08/08/03 (DAK) added drawExpansionFactor to draw in a box larger than
     * simulated one
     * 08/14/03 (DAK) improved selection of default zoom level
     * 04/15/04 (DAK) made init, initialize, and display methods synchronized,
     * to correct problems arising with grand-canonical simulations
     */

//Class used to define canvas onto which configuration is drawn
public class DisplayPhaseCanvas3DOpenGL extends DisplayCanvasOpenGL implements GLUTEnum {
    
  private final double rightClipPlane[] = new double[4];
  private final double leftClipPlane[] = new double[4];
  private final double topClipPlane[] = new double[4];
  private final double bottomClipPlane[] = new double[4];
  private final double frontClipPlane[] = new double[4];
  private final double backClipPlane[] = new double[4];
  private final float MaterialSpecular[] = { 0.8f, 0.8f, 0.8f, 1f };
  private final float MaterialShininess[] = { 70f };
  private final float materialTransparent[] = {0.0f, 0.8f, 0.0f, 0.5f};
  private final float LightSpecular[] = { 1f, 1f, 1f, 1f };
  private final float LightDiffuse[] = { 0.93f, 0.93f, 0.93f, 1f };
  private final float LightPosition[] = { 1f, 1f, 3f, 0f };
  private int sphereList[]; // Storage number for our sphere
  private int wellList[]; // Storage number for our wells
  private int displayList; // Storage number for displaying image shells
  private double sphereListRadius = 0f;
  private double wellListRadius = 0f;
  private boolean glInitialized = false;
  private double drawExpansionFactor = 1.0;
  private final Vector3D vertex = new Vector3D();
  private final Vector3D temp = new Vector3D();

  //The transparent grey color for the wells
  private final static byte wR=(byte)200, wG=(byte)200, wB=(byte)200, wA=(byte)160;
  
  //Variables for translations and zooms
  private float shiftZ = -70f, shiftX = 0f, shiftY = 0f;
  //Rotation Variables
  private float prevx, prevy, xRot = 0f, yRot = 0f;
  //Centers the phase in the canvas
  private float xCenter, yCenter, zCenter;
  private etomica.space3d.Vector3D center = new etomica.space3d.Vector3D();

  private final float vert[];
      
  //Work vector for overflow images
  private float[][] originShifts;
  private double[][] shellOrigins;

  //Local function variables (primarily used in display(), drawDisplay(), and drawBoundary())
  private java.awt.Color lastColor;
  private long T0 = 0, Frames = 0;
  private float drawExpansionShiftX = 0f, drawExpansionShiftY = 0f, drawExpansionShiftZ = 0f;
  //private TextField scaleText = new TextField();
  //private Font font = new Font("sansserif", Font.PLAIN, 10);
  //private int annotationHeight = font.getFontMetrics().getHeight();
  //private int annotationHeight = 12;
      
  public void setShiftX(float x) {shiftX = x;}
  public void setShiftY(float y) {shiftY = y;}
  public void setPrevX(float x) {prevx = x;}
  public void setPrevY(float y) {prevy = y;}
  public void setXRot(float x) {xRot = x;}
  public void setYRot(float y) {yRot = y;}
  public void setZoom(float z) {shiftZ = z;/*System.out.println(z);*/}
  public float getShiftX() {return(shiftX);}
  public float getShiftY() {return(shiftY);}
  public float getPrevX() {return(prevx);}
  public float getPrevY() {return(prevy);}
  public float getXRot() {return(xRot);}
  public float getYRot() {return(yRot);}
  public float getZoom() {return(shiftZ);}
  
  private ColorScheme colorScheme;
  private AtomFilter atomFilter;
          
  public DisplayPhaseCanvas3DOpenGL(DisplayPhase newDisplayPhase, int w, int h) {
    super(w, h);
    //scaleText.setVisible(true);
    //scaleText.setEditable(false);
    //scaleText.setBounds(0,0,100,50);
    vert = new float[3];
    setAnimateFps(24.);
    displayPhase = newDisplayPhase;
    if(displayPhase != null) {
        colorScheme = displayPhase.getColorScheme();
        atomFilter = displayPhase.getAtomFilter();
    }
  }
      
  //public void preInit() {
    //doubleBuffer = true;
    //stereoView = false;
  //}

  public synchronized void init() {
    if(glInitialized) return;
    
    // DAK - 09/27/02 rescale display to adjust to size of phase
    if(displayPhase != null) {
        Phase phase = displayPhase.getPhase();
        if(phase != null) {
            float b = (float)phase.getBoundary().getDimensions().x(0);
//			float z = -70f + (30f/b - 1f)*22f;
//			float z = -190f + (30f/b - 1f)*22f;//08/12/03 DAK changed 70 to 190
			float z = -1.30847f - 2.449f * b;//08/14/03 DAK changed to this by linear regressing observed "best" z vs b values
//			System.out.println(b+"  "+z); 
            setZoom(z);
        }
    }//end DAK
    
    //Set the background clear color
    Color c = getBackground();
    gl.glClearColor((float)c.getRed()/255f, (float)c.getGreen()/255f, (float)c.getBlue()/255f, (float)c.getAlpha()/255f);

    //Enables Clearing Of The Depth Buffer
    gl.glClearDepth(1.0);
    //The Type Of Depth To Do
    gl.glDepthFunc(GL_LEQUAL);
    //Enables Depth Surface Removal
    gl.glEnable(GL_DEPTH_TEST);
    
    //Face culling? worth it or not?
    gl.glEnable(GL_CULL_FACE);
    
    //Disable Dithering, provides speedup on low end systems
    gl.glDisable(GL_DITHER);
    
    //Enable transparency
    gl.glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    gl.glEnable(GL_BLEND);

    //Disable normalization
    gl.glDisable(GL_NORMALIZE);
    
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
 //   	gl.glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//09/17/02


    
    //Let OpenGL know that we wish to use the fastest systems possible
    gl.glHint(GL_CLIP_VOLUME_CLIPPING_HINT_EXT, GL_FASTEST);
    gl.glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    gl.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
    gl.glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    gl.glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
    
    //Set the light properties for the system
    gl.glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
    gl.glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
    gl.glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
    //Enable light
    gl.glEnable(GL_LIGHT0);
    
    //new 09/21/02 (DAK)
    //Set the light properties for the system
 /*   gl.glLightfv(GL_LIGHT1, GL_SPECULAR, LightSpecular);
    gl.glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
    gl.glLightfv(GL_LIGHT1, GL_POSITION, LightPosition2);
    //Enable light
    gl.glEnable(GL_LIGHT1);*/
    //end new
    
    gl.glEnable(GL_LIGHTING);
    
    //Set the material properties for the spheres
    gl.glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, MaterialSpecular);
    gl.glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, MaterialShininess);
    //Set the material to track the current color
    gl.glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    gl.glEnable(GL_COLOR_MATERIAL);

    //Enables Smooth Color Shading
    gl.glShadeModel(GL_SMOOTH);
    
    sphereList = new int[DRAW_QUALITY_MAX];
    wellList = new int[DRAW_QUALITY_MAX];
    sphereList[0] = gl.glGenLists(11);
    wellList[0] = sphereList[0] + 1;
    for(int j = 1; j < DRAW_QUALITY_MAX; j++) {
      sphereList[j] = wellList[j-1] + 1;
      wellList[j] = sphereList[j] + 1;
    }
    displayList = wellList[DRAW_QUALITY_MAX-1] + 1;
    
    setUseRepaint(true);
    setUseFpsSleep(true);
    //setUseRepaint(false);
    //setUseFpsSleep(false);

    rightClipPlane[0] = bottomClipPlane[1] = backClipPlane[2] = 1.;
    leftClipPlane[0] = topClipPlane[1] = frontClipPlane[2] = -1.;
    
	  T0=System.currentTimeMillis();
    start();
    //stop();
    
    glInitialized = true;
  }
  
  public void reshape(int width, int height) {
    if(!glInitialized) return;
    super.reshape(width, height);
    java.awt.Dimension mySize = new java.awt.Dimension(width, height);
	  super.setSize(mySize);
	  setSize(mySize);
    gl.glMatrixMode(GL_PROJECTION);
    gl.glLoadIdentity();
    glu.gluPerspective(35, (float)width/(float)height, 1, 500.0);
    //glu.gluPerspective(45, (float)width/(float)height, 0.1, 100.0);//orig
    gl.glMatrixMode(GL_MODELVIEW);
    gl.glLoadIdentity();
//    gl.glViewport(0,0,width,height);
  }

  private void initSphereList(double sphereRadius) {
    gl.glNewList(sphereList[DRAW_QUALITY_VERY_LOW], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, sphereRadius, 2);
    gl.glEndList();
    gl.glNewList(sphereList[DRAW_QUALITY_LOW], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, sphereRadius, 4);
    gl.glEndList();
    gl.glNewList(sphereList[DRAW_QUALITY_NORMAL], GL_COMPILE);
     gluSphere.gluSmoothSphere(this.gl, sphereRadius, 7);
    gl.glEndList();
    gl.glNewList(sphereList[DRAW_QUALITY_HIGH], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, sphereRadius, 9);
    gl.glEndList();
    gl.glNewList(sphereList[DRAW_QUALITY_VERY_HIGH], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, sphereRadius, 11);
    gl.glEndList();
    sphereListRadius = sphereRadius;
  }
      
  private void initWellList(double wellRadius) {
    gl.glNewList(wellList[DRAW_QUALITY_VERY_LOW], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, wellRadius, 2);
    gl.glEndList();
    gl.glNewList(wellList[DRAW_QUALITY_LOW], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, wellRadius, 4);
    gl.glEndList();
    gl.glNewList(wellList[DRAW_QUALITY_NORMAL], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, wellRadius, 7);
    gl.glEndList();
    gl.glNewList(wellList[DRAW_QUALITY_HIGH], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, wellRadius, 10);
    gl.glEndList();
    gl.glNewList(wellList[DRAW_QUALITY_VERY_HIGH], GL_COMPILE);
    gluSphere.gluSmoothSphere(this.gl, wellRadius, 13);
    gl.glEndList();
    wellListRadius = wellRadius;
  }

  public void setAtomFilter(AtomFilter filter) {atomFilter = filter;}

  private etomica.math.geometry.Plane plane; // (DAK) 09/21/02
  private etomica.space3d.Vector3D[] points;
  private etomica.space3d.Vector3D normal;
  private etomica.space3d.Vector3D nearest;
  protected DisplayPhase displayPhase;
    
  private void drawDisplay() {
    
    colorScheme = displayPhase.getColorScheme();

    if(displayPhase.getColorScheme() instanceof ColorSchemeCollective) {
        ((ColorSchemeCollective)displayPhase.getColorScheme()).colorAllAtoms(displayPhase.getPhase());
    }

    AtomIteratorLeafAtoms atomIter = new AtomIteratorLeafAtoms(displayPhase.getPhase());
    atomIter.reset();

    lastColor = null;
    while (atomIter.hasNext()) {
        AtomLeaf a = (AtomLeaf)atomIter.nextAtom();
        if (a.type instanceof AtomTypeSphere) {
            if(!atomFilter.accept(a)) {continue;}
            
            boolean drawOverflow = false;
            if(displayPhase.getDrawOverflow()) {
                drawOverflow = computeShiftOrigin(a, displayPhase.getPhase().getBoundary());
            }
            Color c = colorScheme.getAtomColor(a);
            Vector r = a.coord.position();
            //Update the positions of the atom
            vert[0] = (float)r.x(0) - xCenter;// + drawExpansionShiftX - 0*(float)temp.x(0);
            vert[1] = (float)r.x(1) - yCenter;// + drawExpansionShiftY - 0*(float)temp.x(1);
            vert[2] = (float)r.x(2) - zCenter;// + drawExpansionShiftZ - 0*(float)temp.x(2);
            //Update the color for the atom
            if(!c.equals(lastColor)) {
              gl.glColor4ub((byte)c.getRed(), (byte)c.getGreen(), (byte)c.getBlue(), (byte)c.getAlpha());
              lastColor = c;
            }
            if(sphereListRadius != ((AtomTypeSphere)a.type).radius(a))
              initSphereList(((AtomTypeSphere)a.type).radius(a));
            gl.glPushMatrix();
            gl.glTranslatef(vert[0], vert[1], vert[2]);
            gl.glCallList(sphereList[getQuality()]);
            //gl.glDrawArrays(GL_TRIANGLE_STRIP, 0, 2);
            //Draw overflow images if so indicated
            if (drawOverflow) {
                int j = originShifts.length;
                while((--j) >= 0) {
                  gl.glPushMatrix();
                  gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
                  gl.glCallList(sphereList[getQuality()]);
                  gl.glPopMatrix();
                }
            }
            gl.glPopMatrix();
            if (a.type instanceof AtomTypeWell) {
                gl.glColor4ub(wR, wG, wB, wA);
                if(wellListRadius != ((AtomTypeWell)a.type).wellRadius())
                    initWellList(((AtomTypeWell)a.type).wellRadius());
                gl.glPushMatrix();
                gl.glTranslatef(vert[0], vert[1], vert[2]);
                gl.glCallList(wellList[getQuality()]);
                //Draw overflow images if so indicated
                if (drawOverflow) {
                    int j = originShifts.length;
                    while((--j) >= 0) {
                      gl.glPushMatrix();
                      gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
                      gl.glCallList(wellList[getQuality()]);
                      gl.glPopMatrix();
                    }
                }
                gl.glPopMatrix();
            }
        }
    }
    
    //Draw other features if indicated
    if(drawBoundary >= DRAW_BOUNDARY_ALL) {
      gl.glDisable(GL_CLIP_PLANE0);
      gl.glDisable(GL_CLIP_PLANE1);
      gl.glDisable(GL_CLIP_PLANE2);
      gl.glDisable(GL_CLIP_PLANE3);
      gl.glDisable(GL_CLIP_PLANE4);
      gl.glDisable(GL_CLIP_PLANE5);
      gl.glDisable(GL_LIGHT0);
      gl.glDisable(GL_LIGHTING);
      gl.glColor3ub((byte)0, (byte)-1, (byte)-1);
      drawBoundary(0);
      gl.glEnable(GL_LIGHT0);
      gl.glEnable(GL_LIGHTING);
//      gl.glEnable(GL_CLIP_PLANE0);
//      gl.glEnable(GL_CLIP_PLANE1);
//      gl.glEnable(GL_CLIP_PLANE2);
//      gl.glEnable(GL_CLIP_PLANE3);
//      gl.glEnable(GL_CLIP_PLANE4);
//      gl.glEnable(GL_CLIP_PLANE5);
    }
    //////////////******drawing of plane*******/////////////////////////
    //(DAK) added this section 09/21/02
    /* do drawing of all drawing objects that have been added to the display */
    for(Iterator iter=displayPhase.getDrawables().iterator(); iter.hasNext(); ) {
      Object obj = iter.next();
      if(obj instanceof etomica.lattice.LatticePlane) {
        plane = ((etomica.lattice.LatticePlane)obj).planeCopy(plane);
        points = plane.inPlaneSquare(plane.nearestPoint(center,nearest),40., points);
        normal = plane.getNormalVector(normal);
    gl.glDisable(GL_CULL_FACE);
        for(int i=0; i<4; i++) points[i].ME(center);
          gl.glBegin(GL_QUADS);
           gl.glNormal3f((float)normal.x(0), (float)normal.x(1), (float)normal.x(2));
//  		   gl.glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
           gl.glMaterialfv(GL_FRONT, GL_AMBIENT, materialTransparent);
           gl.glMaterialfv(GL_FRONT, GL_DIFFUSE, materialTransparent);
           gl.glColor4ub(wR, wG, wB, wA);
           gl.glVertex3f((float)points[0].x(0), (float)points[0].x(1), (float)points[0].x(2));
           gl.glVertex3f((float)points[2].x(0), (float)points[2].x(1), (float)points[2].x(2));
           gl.glVertex3f((float)points[1].x(0), (float)points[1].x(1), (float)points[1].x(2));
           gl.glVertex3f((float)points[3].x(0), (float)points[3].x(1), (float)points[3].x(2));
          gl.glEnd();
        }
    gl.glEnable(GL_CULL_FACE);
    }
    //////////////////////////////////////////
  }
              
  protected boolean computeShiftOrigin(AtomLeaf a, Boundary b) {
    if(a.type instanceof AtomTypeSphere)
      originShifts = b.getOverflowShifts(a.coord.position(),((AtomTypeSphere)a.type).radius(a));
    else
      originShifts = new float[0][0];
    if(originShifts.length == 0) return(false);
    return(true);
  }
        
  /**
  * display is the method that handles the drawing of the phase to the screen.
  */
  public synchronized void display() {
    //Ensure GL is initialized correctly
    if (glj.gljMakeCurrent() == false)
      return;
                            
    //Clear The Screen And The Depth Buffer
    gl.glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    //Reset The View
    gl.glLoadIdentity();
    //PhaseTranslate & Zoom to the desired position
//    gl.glTranslatef(shiftX, shiftY, shiftZ-(displayPhase.getImageShells()<<5));//changed to '5' from '7' (08/12/03 DAK)
	gl.glTranslatef(shiftX, shiftY, (float)(shiftZ-2.5*2*displayPhase.getImageShells()*displayPhase.getPhase().getBoundary().getDimensions().x(0)));//changed to this (08/14/03 DAK)
    //Rotate accordingly
    gl.glRotatef(xRot, 1f, 0f, 0f);
    gl.glRotatef(yRot, 0f, 1f, 0f);
    
    Vector dimensions = displayPhase.getPhase().getBoundary().getDimensions();
    rightClipPlane[3] = leftClipPlane[3] = ((dimensions.x(0))*displayPhase.getImageShells());
    topClipPlane[3] = bottomClipPlane[3] = ((dimensions.x(1))*displayPhase.getImageShells());
    backClipPlane[3] = frontClipPlane[3] = ((dimensions.x(2))*displayPhase.getImageShells());

    gl.glClipPlane(GL_CLIP_PLANE0, rightClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE1, leftClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE2, topClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE3, bottomClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE4, backClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE5, frontClipPlane);

//    gl.glEnable(GL_CLIP_PLANE0);
//    gl.glEnable(GL_CLIP_PLANE1);
//    gl.glEnable(GL_CLIP_PLANE2);
//    gl.glEnable(GL_CLIP_PLANE3);
//    gl.glEnable(GL_CLIP_PLANE4);
//    gl.glEnable(GL_CLIP_PLANE5);

    //Draw periodic images if indicated
    // The following if() block sets up the display list.
    int k = 0;
    if(displayPhase.getImageShells() > 0) {
//      setDrawExpansionFactor(1.0);
      int j = DRAW_QUALITY_VERY_LOW;
      if(displayPhase.getImageShells() == 1) j=DRAW_QUALITY_LOW;
      else if(displayPhase.getImageShells() > 1) j=DRAW_QUALITY_VERY_LOW;
      k = getQuality();
      setQuality(j);
      shellOrigins = displayPhase.getPhase().getBoundary().imageOrigins(displayPhase.getImageShells());
      //more efficient to save rather than recompute each time
      gl.glNewList(displayList, GL_COMPILE_AND_EXECUTE);
    }
    // We always need to draw the display at least once
    drawDisplay();
    // Finish and compile the display list, then redraw it for each shell image
    if(displayPhase.getImageShells() > 0) {
      gl.glEndList();
      int j = shellOrigins.length;
      while((--j) >= 0) {
        gl.glPushMatrix();
        gl.glTranslated(shellOrigins[j][0],shellOrigins[j][1],shellOrigins[j][2]);
        gl.glCallList(displayList);
        gl.glPopMatrix();
      }
      setQuality(k);
    }

    gl.glDisable(GL_CLIP_PLANE0);
    gl.glDisable(GL_CLIP_PLANE1);
    gl.glDisable(GL_CLIP_PLANE2);
    gl.glDisable(GL_CLIP_PLANE3);
    gl.glDisable(GL_CLIP_PLANE4);
    gl.glDisable(GL_CLIP_PLANE5);
    
    //Draw other features if indicated
    if(drawBoundary >= DRAW_BOUNDARY_OUTLINE) {
      gl.glDisable(GL_LIGHT0);
      gl.glDisable(GL_LIGHTING);
      gl.glColor3ub((byte)0, (byte)-1, (byte)-1);
      drawBoundary(displayPhase.getImageShells());
      if (drawBoundary == DRAW_BOUNDARY_SHELL) {
        int j = displayPhase.getImageShells();
        while((--j) >= 0) {
          drawBoundary(j);
        } 
      }
      gl.glEnable(GL_LIGHT0);
      gl.glEnable(GL_LIGHTING);
    }
    

    Frames++;
    long t=System.currentTimeMillis();
    if(t - T0 >= 5000) {
      double seconds = (double)(t - T0) / 1000.0;
      double fps = (double)Frames / seconds;
//      System.out.println(Frames+" frames in "+seconds+" seconds = "+fps+" FPS");
      T0 = t;
      Frames = 0;
    }
    
    //Swap buffers
    glj.gljSwap();
    glj.gljFree();
    //!!!glj.gljCheckGL();
  }
        
  private void drawBoundary(int num) {
    gl.glBegin(GL_LINES);
    Polyhedron shape = (Polyhedron)displayPhase.getPhase().getBoundary().getShape();
    LineSegment[] edges = shape.getEdges();
    for(int i=0; i<edges.length; i++) {
        vertex.E(edges[i].getVertices()[0]);
        vertex.TE((1+2*num));
        gl.glVertex3f((float)vertex.x(0), (float)vertex.x(1), (float)vertex.x(2));
        vertex.E(edges[i].getVertices()[1]);
        vertex.TE((1+2*num));
        gl.glVertex3f((float)vertex.x(0), (float)vertex.x(1), (float)vertex.x(2));
    }
    gl.glEnd();
  }
    
    /**
     * Returns the drawExpansionFactor.
     * @return double
     */
    public double getDrawExpansionFactor() {
        return drawExpansionFactor;
    }
    
    /**
     * Sets the drawExpansionFactor.
     * @param drawExpansionFactor The drawExpansionFactor to set
     */
    public void setDrawExpansionFactor(double drawExpansionFactor) {
        	this.drawExpansionFactor = drawExpansionFactor;
        	if(displayPhase != null && displayPhase.getPhase() != null) {
        		Vector box = displayPhase.getPhase().getBoundary().getDimensions();
        		float mult = (float)(0.5*(drawExpansionFactor - 1.0));//*(2*I+1);
        		drawExpansionShiftX = (float)(mult*box.x(0));
        		drawExpansionShiftY = (float)(mult*box.x(1));
        		drawExpansionShiftZ = (float)(mult*box.x(2));
        	}
    }
}  //end of DisplayPhase.Canvas
