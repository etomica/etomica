//http://www.lanw.com/training/interop/securityurls.htm

//http://www.geocities.com/SiliconValley/Park/5625/opengl/
//http://www.azillionmonkeys.com/windoze/OpenGLvsDirect3D.html

//http://romka.demonews.com/opengl/demos/other_eng.htm
//http://www.oglchallenge.com/

//http://www.sgi.com/software/opengl/advanced97/notes/node196.html#stencilsection
//http://ask.ii.uib.no/ebt-bin/nph-dweb/dynaweb/SGI_Developer/OpenGL_RM
package etomica.graphics;
import etomica.atom.AtomFilter;
import etomica.exception.MethodNotImplementedException;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polyhedron;
import etomica.math.geometry.Polytope;
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
public class DisplayPolytopeCanvas3DOpenGL extends DisplayCanvasOpenGL implements GLUTEnum {
    
  private final double rightClipPlane[] = new double[4];
  private final double leftClipPlane[] = new double[4];
  private final double topClipPlane[] = new double[4];
  private final double bottomClipPlane[] = new double[4];
  private final double frontClipPlane[] = new double[4];
  private final double backClipPlane[] = new double[4];
  private final float LightSpecular[] = { 1f, 1f, 1f, 1f };
  private final float LightDiffuse[] = { 0.93f, 0.93f, 0.93f, 1f };
  private final float LightPosition[] = { 1f, 1f, 3f, 0f };
  private boolean glInitialized = false, canvasInitialized = false;
  private final Vector3D vertex = new Vector3D();

  //Variables for translations and zooms
  private float shiftZ = -70f, shiftX = 0f, shiftY = 0f;
  //Rotation Variables
  private float prevx, prevy, xRot = 0f, yRot = 0f;
  //Centers the phase in the canvas
  private float xCenter, yCenter, zCenter;
  private etomica.space3d.Vector3D center = new etomica.space3d.Vector3D();

  //Local function variables (primarily used in display(), drawDisplay(), and drawBoundary())
//  private long T0 = 0, Frames = 0;
  private java.awt.Color c;
      
  private DisplayPolytope displayPolytope;
  
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
  
  public DisplayPolytopeCanvas3DOpenGL(DisplayPolytope newDisplayPolytope, int w, int h) {
    super(w, h);
    setAnimateFps(24.);
    displayPolytope = newDisplayPolytope;
  }
      
  public void setAtomFilter(AtomFilter foo) {throw new MethodNotImplementedException("I don't wanna talk to you no more, you empty-headed animal food trough water!");}
  
  public synchronized void init() {
    if(glInitialized) return;
    
    // DAK - 09/27/02 rescale display to adjust to size of phase
    if(displayPolytope != null) {
        Polytope polytope = displayPolytope.getPolytope();
        if(polytope != null) {
            float b = (float)displayPolytope.dimensions().x(0);
//			float z = -70f + (30f/b - 1f)*22f;
//			float z = -190f + (30f/b - 1f)*22f;//08/12/03 DAK changed 70 to 190
			float z = -1.30847f - 2.449f * b;//08/14/03 DAK changed to this by linear regressing observed "best" z vs b values
//			System.out.println(b+"  "+z); 
            setZoom(z);
        }
    }//end DAK
    
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
    
    gl.glEnable(GL_LIGHTING);
    
    //Set the material to track the current color
    gl.glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    gl.glEnable(GL_COLOR_MATERIAL);

    //Enables Smooth Color Shading
    gl.glShadeModel(GL_SMOOTH);

    setUseRepaint(true);
    setUseFpsSleep(true);
    //setUseRepaint(false);
    //setUseFpsSleep(false);

    rightClipPlane[0] = bottomClipPlane[1] = backClipPlane[2] = 1.;
    leftClipPlane[0] = topClipPlane[1] = frontClipPlane[2] = -1.;
    
//	  T0=System.currentTimeMillis();
    start();
    //stop();
    
    glInitialized = true;
}
      
  public synchronized void initialize() {
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

  /**
  * Sets the size of the display to a new value and scales the image so that
  * the phase fits in the canvas in the same proportion as before.
  */
  public void scaleSetSize(int width, int height) {
    System.out.println("Scale: New WxH: "+width+"x"+height);
  }
  
  private void drawDisplay() {
    
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
      drawBoundary();
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
              
  /**
  * display is the method that handles the drawing of the phase to the screen.
  */
  public synchronized void display() {
    //Makes sure there is something to draw
    if(!canvasInitialized) {this.initialize();return;}
    //Ensure GL is initialized correctly
    if (glj.gljMakeCurrent() == false)
      return;
                            
    //Clear The Screen And The Depth Buffer
    gl.glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    //Reset The View
    gl.glLoadIdentity();
    //PhaseTranslate & Zoom to the desired position
	gl.glTranslatef(shiftX, shiftY, shiftZ);
    //Rotate accordingly
    gl.glRotatef(xRot, 1f, 0f, 0f);
    gl.glRotatef(yRot, 0f, 1f, 0f);
    
    xCenter = (float)(displayPolytope.dimensions().x(0)*.5);
    yCenter = (float)(displayPolytope.dimensions().x(1)*.5);
    zCenter = (float)(displayPolytope.dimensions().x(2)*.5);
    center.E(xCenter, yCenter, zCenter);
    rightClipPlane[3] = leftClipPlane[3] = xCenter;
    topClipPlane[3] = bottomClipPlane[3] = yCenter;
    backClipPlane[3] = frontClipPlane[3] = zCenter;

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

    // We always need to draw the display at least once
    drawDisplay();

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
      drawBoundary();
      gl.glEnable(GL_LIGHT0);
      gl.glEnable(GL_LIGHTING);
    }
    

//    Frames++;
//    long t=System.currentTimeMillis();
//    if(t - T0 >= 5000) {
//      double seconds = (t - T0) / 1000.0;
//      double fps = Frames / seconds;
//      System.out.println(Frames+" frames in "+seconds+" seconds = "+fps+" FPS");
//      T0 = t;
//      Frames = 0;
//    }
    
    //Swap buffers
    glj.gljSwap();
    glj.gljFree();
    //!!!glj.gljCheckGL();
  }
        
  private void drawBoundary() {
    gl.glBegin(GL_LINES);
    Polyhedron shape = (Polyhedron)displayPolytope.getPolytope();
    LineSegment[] edges = shape.getEdges();
    for(int i=0; i<edges.length; i++) {
        vertex.E(edges[i].getVertices()[0]);
        gl.glVertex3f((float)vertex.x(0), (float)vertex.x(1), (float)vertex.x(2));
        vertex.E(edges[i].getVertices()[1]);
        gl.glVertex3f((float)vertex.x(0), (float)vertex.x(1), (float)vertex.x(2));
    }
    gl.glEnd();
  }

}  //end of DisplayPhase.Canvas
