package g3dsys.images;

import g3dsys.control.G3DSys;
import g3dsys.control.IndexIterator;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

//TODO: implement DRAW_LARGE_BOX type
//endpoints of the lines * (2*numLayers + 1)
//TODO: add lazyBigShell flag to tell us it's okay to scale lines
//canvas sets a flag on the index iterator before passing it
//if flag is true, make large box by scaling boundary lines, else
//hit panic button

public class ImageShell extends Figure {

  //shell boundary drawing modes
  private static final int NUM_DRAW_TYPES = 3; //three types so far
  public static final int DRAW_EVERY_BOX = 0; //redraw boundary for every image
  public static final int DRAW_LARGE_BOX = 1; //draw only one around everything
  public static final int DRAW_INNER_BOX = 2; //draw only the original boundary
  private int drawBoundaryType = 0;
  
  private int numLayers = 1;
  private Graphics3D g3d;
  
  //reusable point objects for conversions
  private Point3f p = new Point3f();
  private Point3f q = new Point3f();
  private Point3i s = new Point3i();
  private Point3i t = new Point3i();
  
  //large box helper
  
  public ImageShell(G3DSys g) {
    super(g,(short)0);
    g3d = g.getG3D();
    drawme = false;
  }
  
  /**
   * Sets the number of image shell layers to use
   * @param l the number of image shell layers
   */
  public void setLayers(int l) { numLayers = (l < 1 ? 1 : l); }
  /**
   * Gets the number of image shell layers
   * @return returns the number of image shell layers
   */
  public int getLayers() { return numLayers; }
  /**
   * Sets the boundary drawing style
   * @param b the boundary drawing style to use
   */
  public void setDrawBoundaryType(int b) {
    //known broken or unimplemented modes should throw exceptions here
    if(b == drawBoundaryType) return; //no change
    if(b >= NUM_DRAW_TYPES) throw new IllegalArgumentException();
    drawBoundaryType = b;
    updateDrawType();
  }
  /**
   * Updates line-drawing status depending on whether images are
   * on and the current draw boundary type. 
   */
  public void updateDrawType() {
    Figure[] figs = _gsys.getFigs();
    for(int i=0; i<figs.length; i++) {
      if(figs[i] == null) continue;
      if(figs[i] instanceof Line) {
        figs[i].setDrawable(!drawme || drawBoundaryType != DRAW_LARGE_BOX);
      }
    }
  }
  /**
   * Gets the boundary drawing style
   * @return returns the boundary drawing style
   */
  public int getDrawBoundaryType() { return drawBoundaryType; } 
  /**
   * Cycles to the next boundary drawing style
   */
  public void cycleDrawBoundaryType() {
    /* Try each successive mode if they fail, up to the maximum of
     * NUM_DRAW_TYPES, then set to zero. This handles the case of one
     * or more broken modes in the middle of the list. Note that this
     * is probably not ever going to happen.
     */ 
    for(int i=1; i<NUM_DRAW_TYPES; i++) {
      try {
        setDrawBoundaryType( (getDrawBoundaryType()+i)%NUM_DRAW_TYPES );
        return;
      } catch(IllegalArgumentException iae) {
        continue;
      }
    }
  }
  
  public void draw() {
    if(!drawme) return;
    double[] vectors = _gsys.getBoundaryVectors();
    int D = vectors.length/3; //dimensionality of boundary; assume 3-dim vectors
    if(D == 0) return; //no vectors, or low dimensionality
    
    //still bound to axis directions; need to use vector orientation instead
    float dx = 0, dy = 0, dz = 0; //for linear combinations
    Figure[] figs = _gsys.getFigs();

    //draw large box if necessary
    if(drawBoundaryType == DRAW_LARGE_BOX) {
      for(int i=0; i<figs.length; i++) {
        if(figs[i] instanceof Line) {
          p.set(((Line)figs[i]).getStart().x * (2*numLayers+1),
              ((Line)figs[i]).getStart().y * (2*numLayers+1),
              ((Line)figs[i]).getStart().z * (2*numLayers+1));
          q.set(((Line)figs[i]).getEnd().x * (2*numLayers+1),
              ((Line)figs[i]).getEnd().y * (2*numLayers+1),
              ((Line)figs[i]).getEnd().z * (2*numLayers+1));
          _gsys.screenSpace(p, s);
          _gsys.screenSpace(q, t);
          g3d.drawDashedLine(figs[i].getColor(),0,0,s.x,s.y,s.z,t.x,t.y,t.z);
        }
      }
    }
    
    //IndexIteratorSequential iter = new IndexIteratorSequential(D,numLayers);
    int[] sizes = new int[D];
    for(int i=0; i<D; i++) { sizes[i] = numLayers*2+1; }
    iter.setSize(sizes);
    iter.reset();
    
    while(iter.hasNext()) {
      int[] ia = iter.next();
      if(ia.length != 3) break; //unexpected dimensionality: bail
      
      //build offsets as linear combination of vectors
      for(int i=0; i<D; i++) {
        dx += (float)( (ia[i]-numLayers) * vectors[3*i]);
        dy += (float)( (ia[i]-numLayers) * vectors[3*i+1]);
        dz += (float)( (ia[i]-numLayers) * vectors[3*i+2]);
      }
      
      if(dx == dy && dy == dz && dz == 0) continue; //don't redraw original
      
      for(int i=0; i<figs.length; i++) {
        Figure f = figs[i];
        if(f == null) continue;
        if(!f.drawme) continue;
        if(f instanceof Ball) {
          if(wireframe) continue; //skip spheres in wireframe mode
          p.set(f.getX()+dx,f.getY()+dy,f.getZ()+dz);
          _gsys.screenSpace(p, s);
          g3d.fillSphereCentered(f.getColor(), _gsys.perspective(s.z, f.getD()), s);
        }
        else if(f instanceof g3dsys.images.Box) {
        }
        else if(f instanceof g3dsys.images.Line) {
          if(drawBoundaryType != DRAW_EVERY_BOX) { continue; }
          p.set(((Line)f).getStart().x+dx,
                ((Line)f).getStart().y+dy,
                ((Line)f).getStart().z+dz);
          q.set(((Line)f).getEnd().x+dx,
                ((Line)f).getEnd().y+dy,
                ((Line)f).getEnd().z+dz);
          _gsys.screenSpace(p, s);
          _gsys.screenSpace(q, t);
          g3d.drawDashedLine(f.getColor(),0,0,s.x,s.y,s.z,t.x,t.y,t.z);
        }
        else if(f instanceof g3dsys.images.Axes) {
        }
        else if(f instanceof g3dsys.images.Bond) {
          Bond b = (Bond) f;
          
          //skip long bonds: assume this means molecule wrapping;
          //may not work as expected for thin systems...
          if( b.getEndpoint1().distance(b.getEndpoint2())
              > _gsys.getAngstromWidth()/2.0f) { continue; }
          
          p.set(b.getEndpoint1().x+dx,
                b.getEndpoint1().y+dy,
                b.getEndpoint1().z+dz);
          q.set(b.getEndpoint2().x+dx,
                b.getEndpoint2().y+dy,
                b.getEndpoint2().z+dz);
          _gsys.screenSpace(p,s);
          _gsys.screenSpace(q,t);
          switch(((Bond)f).getBondType()) {
          case Bond.CYLINDER:
            _gsys.getG3D().fillCylinder(b.getColor1(),b.getColor2(),
                Graphics3D.ENDCAPS_FLAT,(int)b.getD(),
                s.x,s.y,s.z,t.x,t.y,t.z);
            break;
          case Bond.WIREFRAME:
            _gsys.getG3D().drawDashedLine(b.getColor1(), 0, 0, s.x, s.y, s.z, t.x, t.y, t.z);
            //_gsys.getG3D().drawLine(color1,color2,p1i.x,p1i.y,p1i.z,p2i.x,p2i.y,p2i.z);
            break;
          }
        }
      }

      //reset offsets for next iteration
      dx = 0;
      dy = 0;
      dz = 0;
    }
    
  }
  
  public float getD() { return 0; }
  public void setDrawable(boolean b) {
    super.setDrawable(b);
    updateDrawType();
  }

  //local storage of wireframe mode, to propagate effects to images
  private boolean wireframe;
  public void setWireFrame(boolean wireframe) {
    this.wireframe = wireframe;
  }

  private double[] vectors;
  private IndexIterator iter;
  public void setBoundaryVectors(double[] values) { vectors = values; }
  public double[] getBoundaryVectors() { return vectors; }
  public void setBoundaryVectorsIterator(IndexIterator i) { iter = i; }
  public IndexIterator getBoundaryVectorsIterator() { return iter; }
  
}
