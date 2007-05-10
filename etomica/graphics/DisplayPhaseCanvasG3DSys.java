package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Panel;
import java.awt.TextField;

import javax.vecmath.Point3f;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import etomica.space.Boundary;
import etomica.space.IVector;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Bond;
import g3dsys.images.Figure;
import g3dsys.images.Line;

public class DisplayPhaseCanvasG3DSys extends DisplayCanvas
	implements AgentSource, BondManager {

  private TextField scaleText = new TextField();

  //will handle all actual drawing
  private G3DSys gsys;
  private final double[] coords;
	
  private AtomLeafAgentManager aam;
    
  private Polytope oldPolytope;
  private Line[] polytopeLines;
  private boolean boundaryDisplayed = false;
  private Color backgroundColor;
  private Color boundaryFrameColor;
  
  public DisplayPhaseCanvasG3DSys(DisplayPhase _phase) {
    //old stuff
    scaleText.setVisible(true);
    scaleText.setEditable(false);
    scaleText.setBounds(0,0,100,50);
    displayPhase = _phase;

    //init G3DSys
    //adding JPanel flickers, Panel does not. Nobody knows why.
    /*
     * Set visible false here to be toggled later; seems to fix the
     * 'sometimes gray' bug
     */
    //this.setVisible(false); // to be set visible later by SimulationGraphic
    Panel p = new Panel();
    this.setLayout(new java.awt.GridLayout());
    p.setLayout(new java.awt.GridLayout());
    p.setSize(800,800);
    this.add(p);
    coords = new double[3];
    gsys = new G3DSys(p);
    setBackgroundColor(Color.BLACK);
    setBoundaryFrameColor(Color.WHITE);
    //init AtomAgentManager, to sync G3DSys and Etomica models
    //this automatically adds the atoms
    aam = new AtomLeafAgentManager(this, displayPhase.getPhase(), false);

    //gsys.refresh();
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
  
  public void setBackgroundColor(Color color) {
      backgroundColor = color;
      gsys.setBGColor(color);
  }
  
  public Color getBackgroundColor() {
      return backgroundColor;
  }
  
  public void setBoundaryFrameColor(Color color) {
      boundaryFrameColor = color;
      oldPolytope = null;
  }
  
  public Color getBoundaryFrameColor() {
      return boundaryFrameColor;
  }

  public void doPaint(Graphics g) {
    
    //handle pending bond addition requests
    if(pendingBonds.size() > 0) {
      for(int i=0; i<pendingBonds.size(); i++) {
        Object[] o = (Object[]) pendingBonds.get(i);
        Ball ball0 = (Ball) o[0];
        Ball ball1 = (Ball) o[1];
        //can't do anything with bondType for now
        Figure f = new Bond(gsys,ball0,ball1);
        gsys.addFig(f);
      }
    }
    
    ColorScheme colorScheme = displayPhase.getColorScheme();
    AtomFilter atomFilter = displayPhase.getAtomFilter();
    if(colorScheme instanceof ColorSchemeCollective) {
      ((ColorSchemeCollective)colorScheme).colorAllAtoms();
    }
    AtomArrayList leafList = displayPhase.getPhase().getSpeciesMaster().getLeafList();
    int nLeaf = leafList.size();
    for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
        IAtomPositioned a = (IAtomPositioned)leafList.get(iLeaf);
      if (a==null || !(a.getType() instanceof AtomTypeSphere)) continue;
      Ball ball = (Ball)aam.getAgent(a);
      if (ball == null) {
        continue;
      }
      /*
       * Atomfilter changes the drawable flag in spheres; bonds respect
       * this and will not draw themselves either. Wireframe mode, on the
       * other hand, tells G3DSys to ignore spheres entirely regardless
       * of drawable flag. This makes it possible to filter bonds in
       * wireframe mode as well.
       */
      boolean drawable = atomFilter.accept(a);
      ball.setDrawable(drawable);
      if (!drawable) {
        continue;
      }
      a.getPosition().assignTo(coords);
      float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
      ball.setColor(G3DSys.getColix(colorScheme.getAtomColor(a)));
      ball.setD(diameter);
      ball.setX((float)coords[0]);
      ball.setY((float)coords[1]);
      ball.setZ((float)coords[2]);
    }
        
    Boundary boundary = displayPhase.getPhase().getBoundary();
    Polytope polytope = boundary.getShape();
    if (polytope != oldPolytope) {

      //force trunc. oct. to make vecs else null pointer exception
      boundary.getPeriodicVectors(); 
      //send iterator to g3dsys
      gsys.setBoundaryVectorsIterator(
          wrapIndexIterator((boundary.getIndexIterator())));

      if (polytopeLines != null) {
        for (int i=0; i<polytopeLines.length; i++) {
          gsys.removeFig(polytopeLines[i]);
        }
      }
      LineSegment[] lines = polytope.getEdges();
      polytopeLines = new Line[lines.length];
      for (int i=0; i<lines.length; i++) {
        IVector[] vertices = lines[i].getVertices();
        polytopeLines[i] = new Line(gsys, G3DSys.getColix(boundaryFrameColor), 
            new Point3f((float)vertices[0].x(0), (float)vertices[0].x(1), (float)vertices[0].x(2)), 
            new Point3f((float)vertices[1].x(0), (float)vertices[1].x(1), (float)vertices[1].x(2)));
        if (displayPhase.getShowBoundary() == true) {
            gsys.addFig(polytopeLines[i]);
        }
      }
      oldPolytope = polytope;
    }
    else {
      LineSegment[] lines = polytope.getEdges();
      for (int i=0; i<lines.length; i++) {
        IVector[] vertices = lines[i].getVertices();
        polytopeLines[i].setStart((float)vertices[0].x(0), (float)vertices[0].x(1), (float)vertices[0].x(2));
        polytopeLines[i].setEnd((float)vertices[1].x(0), (float)vertices[1].x(1), (float)vertices[1].x(2));

        if (displayPhase.getShowBoundary() == false &&
        	boundaryDisplayed == true) {
        	gsys.removeFig(polytopeLines[i]);
        }
        else if (displayPhase.getShowBoundary() == true &&
        		 boundaryDisplayed == false) {
        	gsys.addFig(polytopeLines[i]);
        }
      }
    }

    if (displayPhase.getShowBoundary() == false) {
  	  boundaryDisplayed = false;
    }
    else {
  	  boundaryDisplayed = true;
    }

    // set boundary vectors for image shell
    IVector[] vecs = boundary.getPeriodicVectors();
    double[] dvecs = new double[vecs.length*3]; //assuming 3-dimensional vectors
    for(int i=0; i<vecs.length; i++) {
      if(vecs[i] == null) continue;
      dvecs[i*3] = vecs[i].x(0);
      dvecs[i*3+1] = vecs[i].x(1);
      dvecs[i*3+2] = vecs[i].x(2);
    }
    gsys.setBoundaryVectors(dvecs);
    
    IVector bounds = boundary.getBoundingBox();
    gsys.setBoundingBox((float)(-bounds.x(0)*0.5), (float)(-bounds.x(1)*0.5), (float)(-bounds.x(2)*0.5),
        (float)( bounds.x(0)*0.5), (float)( bounds.x(1)*0.5), (float)( bounds.x(2)*0.5));
          
    gsys.fastRefresh();
  }
  
  /**
   * Add a bond to the graphical display between the given pairs.  The given
   * bondType is used to decide how the bond should be drawn.
   */
  public Object makeBond(AtomSet pair, Object bondType) {
    /*
     * Ball objects here could be null if the bond is created before
     * the atoms have been added. Check for this and store atoms locally
     * in a list. In doPaint check list for pending additions and add
     * them.
     */
    //bondType is a potential right now
    //best to ignore it for now; all bonds are equal
    Ball ball0 = (Ball)aam.getAgent(pair.getAtom(0));
    Ball ball1 = (Ball)aam.getAgent(pair.getAtom(1));
    if(ball0 == null || ball1 == null) {
      System.out.println("NULL!!!");
      pendingBonds.add(new Object[] { ball0, ball1, bondType } );
      return null;
    }
    
    // make a bond object (Figure)
    Figure f = new Bond(gsys,ball0,ball1);
    gsys.addFig(f);
    return f;
  }
  
  private java.util.ArrayList pendingBonds = new java.util.ArrayList();
  
  /**
   * Removes the given bond from the graphical display.  The bond must be an
   * Object returned by the makeBond method.
   */
  public void releaseBond(Object bond) {
    Figure figure = (Figure)bond;
    if (figure.getID() == -1) {
      throw new RuntimeException(figure+" has already been removed");
    }
    gsys.removeFig(figure);
  }
    
  /* ******************************************************
   * AgentSource methods
   * ******************************************************/
  public Class getAgentClass() {
    return Figure.class;
  }
	
  public Object makeAgent(IAtom a) {
    if (!(a.getType() instanceof AtomTypeSphere)) return null;
    ((IAtomPositioned)a).getPosition().assignTo(coords);
    
    float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
    Ball newBall = new Ball(gsys, G3DSys.getColix((displayPhase.getColorScheme().getAtomColor(a))),
        (float)coords[0], (float)coords[1], (float)coords[2], diameter);
    gsys.addFig(newBall);
    return newBall;
  }

  public void releaseAgent(Object agent, IAtom atom) {
    gsys.removeFig((Figure) agent);
  }
  
  /**
   * Set slab percentage
   * @param slab the slab percentage to set
   */
  public void setSlab(double slab) { gsys.setSlabPercent((int)slab); }
  /**
   * Get depth percentage
   * @return returns current depth percentage
   */
  public double getSlab() { return gsys.getSlabPercent(); }
  /**
   * Set depth percentage
   * @param depth the depth percentage to set
   */
  public void setDepth(double depth) { gsys.setDepthPercent((int)depth); }
  /**
   * Get depth percentage
   * @return returns current depth percentage
   */
  public double getDepth() { return gsys.getDepthPercent(); }
  

  /**
   * Wraps an etomica index iterator in an equivalent g3dsys interface
   * for transport; removes g3dsys dependency from all but the
   * etomica.graphics package.
   * @param iter the etomica index iterator to wrap
   * @return returns the g3dsys index iterator
   */
  private g3dsys.control.IndexIterator wrapIndexIterator(
      etomica.lattice.IndexIteratorSizable iter) {
    
    final etomica.lattice.IndexIteratorSizable i = iter;
    
    return new g3dsys.control.IndexIterator() {

      private etomica.lattice.IndexIteratorSizable ii = i;
      
      public int getD() { return ii.getD(); }
      public boolean hasNext() { return ii.hasNext(); }
      public int[] next() { return ii.next(); }
      public void reset() { ii.reset(); }
      public void setSize(int[] size) { ii.setSize(size); }
      
      public boolean isLazySafe() {
        /*
         * For now all boundaries are lazy-safe, including truncated
         * octahedron. If this changes, check for that boundary type
         * here (instanceof IndexIteratorSequentialFiltered, say,
         * after making the class public) and use appropriate boolean.
         */
        return true;
      }
    
    };
    
  }

  
}
