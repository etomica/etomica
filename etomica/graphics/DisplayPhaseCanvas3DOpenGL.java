//http://www.lanw.com/training/interop/securityurls.htm

//http://www.geocities.com/SiliconValley/Park/5625/opengl/
//http://www.azillionmonkeys.com/windoze/OpenGLvsDirect3D.html

//http://romka.demonews.com/opengl/demos/other_eng.htm
//http://www.oglchallenge.com/

//http://www.sgi.com/software/opengl/advanced97/notes/node196.html#stencilsection
//http://ask.ii.uib.no/ebt-bin/nph-dweb/dynaweb/SGI_Developer/OpenGL_RM
package etomica.graphics;
import etomica.*;

import gl4java.utils.glut.fonts.data.*;
import gl4java.utils.glut.fonts.*;
import gl4java.utils.glut.*;

import java.awt.Font;
import java.awt.Shape;
import java.awt.font.*;
import java.awt.geom.*;

import etomica.utility.Iterator;

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
  private final float LightSpecular[] = { 1f, 1f, 1f, 1f };
  private final float LightDiffuse[] = { 0.93f, 0.93f, 0.93f, 1f };
  private final float LightPosition[] = { 1f, 1f, 3f, 0f };
  private int sphereList[]; // Storage number for our sphere
  private int wellList[]; // Storage number for our wells
  private int displayList; // Storage number for displaying image shells
  private double sphereListRadius = 0f;
  private double wellListRadius = 0f;
  private boolean glInitialized = false, canvasInitialized = false;

  //The transparent grey color for the wells
  private final static byte wR=(byte)200, wG=(byte)200, wB=(byte)200, wA=(byte)160;
  //The marker color for the orientations
  private final static byte mR=(byte)255, mG=(byte)10, mB=(byte)60;
  
  //Variables for translations and zooms
  private float shiftZ = -70f, shiftX = 0f, shiftY = 0f;
  //Rotation Variables
  private float prevx, prevy, xRot = 0f, yRot = 0f;
  //Centers the phase in the canvas
  private float xCenter, yCenter, zCenter;

  //The groups of atoms
  private Atom sphereCores[];
  private Atom sphereWells[];
  private Atom sphereRotators[];
  private Atom walls[];
  //The verticies of said atoms
  private float vertSphereCores[];
  private int vertSphereWellBase[];
  private int vertSphereRotatorBase[];
  private float vertWalls[];
      
  //Work vector for overflow images
  private float[][] originShifts;
  private double[][] shellOrigins;

  //Local function variables (primarily used in display(), drawDisplay(), and drawBoundary())
  private float xCent, yCent, zCent;
  private java.awt.Color lastColor;
  private long T0 = 0, Frames = 0;
  private java.awt.Color c;
  private etomica.Atom a;
  private Space.Vector r;
  private int i, j, k;
  
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
  public void setZoom(float z) {shiftZ = z;}
  public float getShiftX() {return(shiftX);}
  public float getShiftY() {return(shiftY);}
  public float getPrevX() {return(prevx);}
  public float getPrevY() {return(prevy);}
  public float getXRot() {return(xRot);}
  public float getYRot() {return(yRot);}
  public float getZoom() {return(shiftZ);}
  
  private ColorScheme colorScheme;
          
  public DisplayPhaseCanvas3DOpenGL(DisplayPhase newDisplayPhase, int w, int h) {
    super(w, h);
    //scaleText.setVisible(true);
    //scaleText.setEditable(false);
    //scaleText.setBounds(0,0,100,50);
    setAnimateFps(24.);
    displayPhase = newDisplayPhase;
    if(displayPhase != null) colorScheme = displayPhase.getColorScheme();
  }
      
  //public void preInit() {
    //doubleBuffer = true;
    //stereoView = false;
  //}

  public void init() {
    if(glInitialized) return;
    //Set the background clear color
    c = getBackground();
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
    for(j = 1; j < DRAW_QUALITY_MAX; j++) {
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
      
  public void initialize() {
    int countSphereCores = 0, countSphereWells = 0, countSphereRotators = 0;
    int countWalls = 0, countAll = 0;
    float vertAll[];
    Atom atoms[];

    countAll = displayPhase.getPhase().speciesMaster.leafAtomCount();
    
    if(countAll==0) return;
    
    vertAll = new float[countAll*3];
    atoms = new Atom[countAll];
    
    i = 0;
    for(Atom a = displayPhase.getPhase().firstAtom(); a!=null; a=a.nextAtom(), i+=3) {
      atoms[i/3] = a;
      vertAll[i] = (float)a.coord.position().component(0);
      vertAll[i+1] = (float)a.coord.position().component(1);
      vertAll[i+2] = (float)a.coord.position().component(2);
      if(a.type instanceof AtomType.OrientedSphere) countSphereRotators++;
      if(a.type instanceof AtomType.Well) countSphereWells++;
      if(a.type instanceof AtomType.Sphere) countSphereCores++;
      if(a.type instanceof AtomType.Wall) countWalls++;
    }
    
    sphereCores = new Atom[countSphereCores];
    sphereWells = new Atom[countSphereWells];
    sphereRotators = new Atom[countSphereRotators];
    walls = new Atom[countWalls];
    vertSphereCores = new float[countSphereCores*3];
    vertSphereWellBase = new int[countSphereWells];
    vertSphereRotatorBase = new int[countSphereRotators];
    vertWalls = new float[countWalls*3];
    for(int j=0,k=0,l=0,m=0,n=0; (j/3) < atoms.length; j+=3) {
      if(atoms[j/3].type instanceof AtomType.Sphere) {
        sphereCores[m/3] = atoms[j/3];
        vertSphereCores[m] = vertAll[j];
        vertSphereCores[m+1] = vertAll[j+1];
        vertSphereCores[m+2] = vertAll[j+2];
        m+=3;
      }
      if(atoms[j/3].type instanceof AtomType.OrientedSphere) {
        sphereRotators[k] = atoms[j/3];
        vertSphereRotatorBase[k] = m-3;
        k++;
      }
      if(atoms[j/3].type instanceof AtomType.Well) {
        sphereWells[l] = atoms[j/3];
        vertSphereWellBase[l] = m-3;
        l++;
      }
      if(atoms[j/3].type instanceof AtomType.Wall) {
        walls[n/3] = atoms[j/3];
        vertWalls[n] = vertAll[j];
        vertWalls[n+1] = vertAll[j+1];
        vertWalls[n+2] = vertAll[j+2];
        n+=3;
      }
    }
    canvasInitialized = true;
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

  /**
  * Sets the size of the display to a new value and scales the image so that
  * the phase fits in the canvas in the same proportion as before.
  */
  public void scaleSetSize(int width, int height) {
    System.out.println("Scale: New WxH: "+width+"x"+height);
    /*if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
      double ratio1 = (double)width/(double)getBounds().width;
      double ratio2 = (double)height/(double)getBounds().height;
      double factor = Math.min(ratio1, ratio2);
      //double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
      displayPhase.setScale(displayPhase.getScale()*factor);
      setSize(width, height);
    }*/
  }

  private void drawDisplay() {

    /* do drawing of all drawing objects that have been added to the display */
    /*for(java.util.Iterator iter=displayPhase.getDrawables().iterator(); iter.hasNext(); ) {
      Drawable obj = (Drawable)iter.next();
    }*/

    if(walls.length > 0) {
      lastColor = null;
      i = walls.length - 1;
      i += i<<1;
      while(i >= 0) {
        a = walls[i/3];
        c = colorScheme.atomColor(a);
        r = a.coord.position();
        //Update the positions of the atom
        vertWalls[i] = (float)r.component(0) - xCenter;
        vertWalls[i+1] = (float)r.component(1) - yCenter;
        vertWalls[i+2] = (float)r.component(2) - zCenter;
        //Update the color for the atom
        if(!c.equals(lastColor)) {
          gl.glColor4ub((byte)c.getRed(), (byte)c.getGreen(), (byte)c.getBlue(), (byte)c.getAlpha());
          lastColor = c;
        }
        gl.glPushMatrix();
        gl.glTranslatef(vertWalls[i], vertWalls[i+1], vertWalls[i+2]);
        //!!! Draw wall here
        //Draw overflow images if so indicated
        if(displayPhase.getDrawOverflow()) {
          if(computeShiftOrigin(a, displayPhase.getPhase().boundary())) {
            j = originShifts.length;
            while((--j) >= 0) {
              gl.glPushMatrix();
              gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
              //!!! Draw wall here
              gl.glPopMatrix();
            }
          }
        }
        gl.glPopMatrix();
        i-=3;
      }
    }
    if(sphereCores.length > 0) {
      lastColor = null;
      i = sphereCores.length - 1;
      i += i<<1;
      while(i >= 0) {
        a = sphereCores[i/3];
        Space3D.Vector r3 = (Space3D.Vector)a.coord.position();
//        System.out.println(a.toString()+", "+r3.x+", "+r3.y+", "+r3.z);
        c = colorScheme.atomColor(a);
        r = a.coord.position();
        //Update the positions of the atom
        vertSphereCores[i] = (float)r.component(0) - xCenter;
        vertSphereCores[i+1] = (float)r.component(1) - yCenter;
        vertSphereCores[i+2] = (float)r.component(2) - zCenter;
        //Update the color for the atom
        if(!c.equals(lastColor)) {
          gl.glColor4ub((byte)c.getRed(), (byte)c.getGreen(), (byte)c.getBlue(), (byte)c.getAlpha());
          lastColor = c;
        }
        if(sphereListRadius != ((AtomType.Sphere)a.type).radius())
          initSphereList(((AtomType.Sphere)a.type).radius());
        gl.glPushMatrix();
        gl.glTranslatef(vertSphereCores[i], vertSphereCores[i+1], vertSphereCores[i+2]);
        gl.glCallList(sphereList[getQuality()]);
        //gl.glDrawArrays(GL_TRIANGLE_STRIP, 0, 2);
        //Draw overflow images if so indicated
        if(displayPhase.getDrawOverflow()) {
          if(computeShiftOrigin(a, displayPhase.getPhase().boundary())) {
            j = originShifts.length;
            while((--j) >= 0) {
              gl.glPushMatrix();
              gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
              gl.glCallList(sphereList[getQuality()]);
              gl.glPopMatrix();
            }
          }
        }
        gl.glPopMatrix();
        i-=3;
      }
    }
    if(sphereWells.length > 0) {
      i = sphereWells.length;
      gl.glColor4ub(wR, wG, wB, wA);
      while((--i) >= 0) {
        a = sphereWells[i];
        if(wellListRadius != ((AtomType.Well)a.type).wellRadius())
          initWellList(((AtomType.Well)a.type).wellRadius());
        gl.glPushMatrix();
        gl.glTranslatef(vertSphereCores[vertSphereWellBase[i]], vertSphereCores[vertSphereWellBase[i]+1], vertSphereCores[vertSphereWellBase[i]+2]);
        gl.glCallList(wellList[getQuality()]);
        //Draw overflow images if so indicated
        if(displayPhase.getDrawOverflow()) {
          if(computeShiftOrigin(a, displayPhase.getPhase().boundary())) {
            j = originShifts.length;
            while((--j) >= 0) {
              gl.glPushMatrix();
              gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
              gl.glCallList(wellList[getQuality()]);
              gl.glPopMatrix();
            }
          }
        }
        gl.glPopMatrix();
      }
    }
    if(sphereRotators.length > 0) {
      i = sphereRotators.length;
      gl.glColor3ub(mR, mG, mB);
      while((--i) >= 0) {
        a = sphereRotators[i];
        gl.glPushMatrix();
        gl.glTranslatef(vertSphereCores[vertSphereRotatorBase[i]], vertSphereCores[vertSphereRotatorBase[i]+1], vertSphereCores[vertSphereRotatorBase[i]+2]);
        ///!!! Draw rotator orientation here
        //Draw overflow images if so indicated
        if(displayPhase.getDrawOverflow()) {
          if(computeShiftOrigin(a, displayPhase.getPhase().boundary())) {
            j = originShifts.length;
            while((--j) >= 0) {
              gl.glPushMatrix();
              gl.glTranslatef(originShifts[j][0], originShifts[j][1], originShifts[j][2]);
              ///!!! Draw rotator orientation here
              gl.glPopMatrix();
            }
          }
        }
        gl.glPopMatrix();
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
      gl.glEnable(GL_CLIP_PLANE0);
      gl.glEnable(GL_CLIP_PLANE1);
      gl.glEnable(GL_CLIP_PLANE2);
      gl.glEnable(GL_CLIP_PLANE3);
      gl.glEnable(GL_CLIP_PLANE4);
      gl.glEnable(GL_CLIP_PLANE5);
    }
  }
              
  protected boolean computeShiftOrigin(Atom a, Space.Boundary b) {
    if(a.type instanceof AtomType.Sphere)
      originShifts = b.getOverflowShifts(a.coord.position(),((AtomType.Sphere)a.type).radius());
    else if(a.type instanceof AtomType.Wall)
      originShifts = b.getOverflowShifts(a.coord.position(),((AtomType.Wall)a.type).getLength());
    else
      originShifts = new float[0][0];
    if(originShifts.length == 0) return(false);
    return(true);
  }
        
  /**
  * display is the method that handles the drawing of the phase to the screen.
  */
  public void display() {
    //Makes sure there is something to draw
    if(!canvasInitialized) {this.initialize();return;}
    //Ensure GL is initialized correctly
    if (glj.gljMakeCurrent() == false)
      return;
                            
    //Clear The Screen And The Depth Buffer
    gl.glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    //Reset The View
    gl.glLoadIdentity();
    //Translate & Zoom to the desired position
    gl.glTranslatef(shiftX, shiftY, shiftZ-(displayPhase.getImageShells()<<7));
    //Rotate accordingly
    gl.glRotatef(xRot, 1f, 0f, 0f);
    gl.glRotatef(yRot, 0f, 1f, 0f);
    
    //Color all atoms according to colorScheme in DisplayPhase
//    displayPhase.getColorScheme().colorAllAtoms();

    xCenter = (float)(displayPhase.getPhase().boundary().dimensions().component(0)*.5);
    yCenter = (float)(displayPhase.getPhase().boundary().dimensions().component(1)*.5);
    zCenter = (float)(displayPhase.getPhase().boundary().dimensions().component(2)*.5);
    rightClipPlane[3] = leftClipPlane[3] = xCenter + ((2*xCenter)*displayPhase.getImageShells());
    topClipPlane[3] = bottomClipPlane[3] = yCenter + ((2*yCenter)*displayPhase.getImageShells());
    backClipPlane[3] = frontClipPlane[3] = zCenter + ((2*zCenter)*displayPhase.getImageShells());

    gl.glClipPlane(GL_CLIP_PLANE0, rightClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE1, leftClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE2, topClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE3, bottomClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE4, backClipPlane);
    gl.glClipPlane(GL_CLIP_PLANE5, frontClipPlane);

    gl.glEnable(GL_CLIP_PLANE0);
    gl.glEnable(GL_CLIP_PLANE1);
    gl.glEnable(GL_CLIP_PLANE2);
    gl.glEnable(GL_CLIP_PLANE3);
    gl.glEnable(GL_CLIP_PLANE4);
    gl.glEnable(GL_CLIP_PLANE5);

    //Draw periodic images if indicated
    // The following if() block sets up the display list.
    if(displayPhase.getImageShells() > 0) {
      if(displayPhase.getImageShells() == 1) j=DRAW_QUALITY_LOW;
      else if(displayPhase.getImageShells() > 1) j=DRAW_QUALITY_VERY_LOW;
      k = getQuality();
      setQuality(j);
      shellOrigins = displayPhase.getPhase().boundary().imageOrigins(displayPhase.getImageShells());
      //more efficient to save rather than recompute each time
      gl.glNewList(displayList, GL_COMPILE_AND_EXECUTE);
    }
    // We always need to draw the display at least once
    drawDisplay();
    // Finish and compile the display list, then redraw it for each shell image
    if(displayPhase.getImageShells() > 0) {
      gl.glEndList();
      j = shellOrigins.length;
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
        j = displayPhase.getImageShells();
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
      System.out.println(Frames+" frames in "+seconds+" seconds = "+fps+" FPS");
      T0 = t;
      Frames = 0;
    }
    
    //Swap buffers
    glj.gljSwap();
    glj.gljFree();
    //!!!glj.gljCheckGL();
  }
        
  private void drawBoundary(int num) {
    xCent = xCenter+((2*xCenter)*num);
    yCent = yCenter+((2*yCenter)*num);
    zCent = zCenter+((2*zCenter)*num);
    gl.glBegin(GL_LINES);
    gl.glVertex3f(-xCent, yCent, zCent);
    gl.glVertex3f(xCent, yCent, zCent);
    gl.glVertex3f(xCent, yCent, zCent);
    gl.glVertex3f(xCent, -yCent, zCent);
    gl.glVertex3f(xCent, -yCent, zCent);
    gl.glVertex3f(-xCent, -yCent, zCent);
    gl.glVertex3f(-xCent, -yCent, zCent);
    gl.glVertex3f(-xCent, -yCent, -zCent);
    gl.glVertex3f(-xCent, -yCent, -zCent);
    gl.glVertex3f(xCent, -yCent, -zCent);
    gl.glVertex3f(xCent, -yCent, -zCent);
    gl.glVertex3f(xCent, yCent, -zCent);
    gl.glVertex3f(xCent, yCent, -zCent);
    gl.glVertex3f(-xCent, yCent, -zCent);
    gl.glVertex3f(-xCent, yCent, -zCent);
    gl.glVertex3f(-xCent, yCent, zCent);
    gl.glVertex3f(xCent, yCent, zCent);
    gl.glVertex3f(xCent, yCent, -zCent);
    gl.glVertex3f(-xCent, -yCent, -zCent);
    gl.glVertex3f(-xCent, yCent, -zCent);
    gl.glVertex3f(xCent, -yCent, zCent);
    gl.glVertex3f(xCent, -yCent, -zCent);
    gl.glVertex3f(-xCent, yCent, zCent);
    gl.glVertex3f(-xCent, -yCent, zCent);
    gl.glEnd();
  }
    
  /*public void doPaint(Graphics g) {
    //Draw bar showing scale if indicated
    if(writeScale) {
      g.setColor(Color.lightGray);
      g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
      g.setColor(Color.black);
      g.setFont(font);
      g.drawString("Scale: "+Integer.toString((int)(100*displayPhase.getScale()))+"%", 0, getSize().height-3);
    }
    //System.gc();
  }*/
}  //end of DisplayPhase.Canvas
