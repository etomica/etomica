package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

import g3dsys.control.G3DSys;

public class ImageShell extends Figure {

  private int numLayers = 4;
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
    Figure[] figs = _gsys.getFigs();
    float dx = _gsys.getAngstromWidth();
    float dy = _gsys.getAngstromHeight();
    float dz = _gsys.getAngstromHeight();
    
    for(int z=-numLayers; z<=numLayers; z++) { //front to back...
      for(int y=-numLayers; y<=numLayers; y++) { //top to bottom...
        for(int x=-numLayers; x<=numLayers; x++) { //left to right...
          if(x == y && y == z && z == 0) continue; //don't redraw original
          
          for(int i=0; i<figs.length; i++) {
            Figure f = figs[i];
            if(f == null) continue;
            if(!f.drawme) continue;
            if(f instanceof Ball) {
              p.set(f.getX()+(dx*x),f.getY()+(dy*y),f.getZ()+(dz*z));
              _gsys.screenSpace(p, s);
              g3d.fillSphereCentered(f.getColor(), _gsys.perspective(s.z, f.getD()), s);
            }
            else if(f instanceof g3dsys.images.Box) {
            }
            else if(f instanceof g3dsys.images.Line) {
              p.set(((Line)f).getStart().x+(dx*x),
                    ((Line)f).getStart().y+(dy*y),
                    ((Line)f).getStart().z+(dz*z));
              q.set(((Line)f).getEnd().x+(dx*x),
                    ((Line)f).getEnd().y+(dy*y),
                    ((Line)f).getEnd().z+(dz*z));
              _gsys.screenSpace(p, s);
              _gsys.screenSpace(q, t);
              g3d.drawDashedLine(f.getColor(),0,0,s.x,s.y,s.z,t.x,t.y,t.z);
            }
            else if(f instanceof g3dsys.images.Axes) {
            }
            else if(f instanceof g3dsys.images.Cylinder) {
            }
          }
          
        }
      }
    }
    
  }

  public float getD() { return 0; }

}
