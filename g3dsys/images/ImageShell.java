package g3dsys.images;

import g3dsys.control.G3DSys;
import g3dsys.control.IndexIteratorSequential;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

public class ImageShell extends Figure {

  private int numLayers = 1;
  private Graphics3D g3d;
  
  //reusable point objects for conversions
  private Point3f p = new Point3f();
  private Point3f q = new Point3f();
  private Point3i s = new Point3i();
  private Point3i t = new Point3i();
  
  public ImageShell(G3DSys g) {
    super(g,(short)0);
    g3d = g.getG3D(); 
  }
  
  public void setLayers(int l) { numLayers = l; }
  public int getLayers() { return numLayers; }
  
  public void draw() {
    double[] vectors = _gsys.getBoundaryVectors();
    int D = vectors.length/3; //dimensionality of boundary; assume 3-dim vectors
    if(D == 0) return; //no vectors, or low dimensionality
    
    //still bound to axis directions; need to use vector orientation instead
    float dx = 0, dy = 0, dz = 0; //for linear combinations
    Figure[] figs = _gsys.getFigs();
    
    IndexIteratorSequential iter = new IndexIteratorSequential(D,numLayers);
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
          p.set(f.getX()+dx,f.getY()+dy,f.getZ()+dz);
          _gsys.screenSpace(p, s);
          g3d.fillSphereCentered(f.getColor(), _gsys.perspective(s.z, f.getD()), s);
        }
        else if(f instanceof g3dsys.images.Box) {
        }
        else if(f instanceof g3dsys.images.Line) {
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
        else if(f instanceof g3dsys.images.Cylinder) {
        }
      }

      //reset offsets for next iteration
      dx = 0;
      dy = 0;
      dz = 0;
    }
    
  }
  
  public float getD() { return 0; }

}
