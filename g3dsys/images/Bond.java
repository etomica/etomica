package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

import g3dsys.control.G3DSys;

public class Bond extends Figure {

  public static final int CYLINDER = 0;
  public static final int WIREFRAME = 1;
  
  private final Point3f p1; // center of first endcap
  private final Point3f p2; // center of second endcap
  private final Point3i p1i, p2i; // same in pixels
  private int bondType = CYLINDER;
  private short color1;
  private short color2;
  
  public Bond(G3DSys g, short c, float x1, float y1, float z1, float x2,
      float y2, float z2) {
    this(g, c, c, new Point3f(x1, y1, z1), new Point3f(x2, y2, z2));
  }
  
  public Bond(G3DSys g, short c1, short c2, float x1, float y1, float z1,
      float x2, float y2, float z2) {
    this(g, c1, c2, new Point3f(x1, y1, z1), new Point3f(x2, y2, z2));
  }


  public Bond(G3DSys g, short c, short c2, Point3f p1, Point3f p2) {
    super(g, c);
    this.p1 = p1;
    this.p2 = p2;
    _p = new Point3f((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
    p1i = new Point3i();
    p2i = new Point3i();
    color1 = c;
    color2 = c2;
  }

  public void draw() {
    _gsys.screenSpace(p1, p1i);
    _gsys.screenSpace(p2, p2i);
    switch(bondType) {
    case CYLINDER:
      _gsys.getG3D().fillCylinder(color1,color2,Graphics3D.ENDCAPS_FLAT,
          (int)getD(),p1i.x,p1i.y,p1i.z,p2i.x,p2i.y,p2i.z);
      break;
    case WIREFRAME:
      _gsys.getG3D().drawDashedLine(color1, 0, 0, p1i.x, p1i.y, p1i.z, p2i.x, p2i.y, p2i.z);
      //_gsys.getG3D().drawLine(color1,color2,p1i.x,p1i.y,p1i.z,p2i.x,p2i.y,p2i.z);
      break;
    }

  }

  public float getD() {
    // no good way to pick diameter yet
    // also, what does this mean when wireframe is on?
    return p1.distance(p2) / 3.0f * (_gsys.getPixelWidth() / _gsys.getAngstromWidth());
  }

}
