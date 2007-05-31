package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;

import g3dsys.control.G3DSys;

//TODO: ensure wireframe bonds actually connect

public class Bond extends Figure {

  public static final int CYLINDER = 0;
  public static final int WIREFRAME = 1;
  
  //points are shared with the Ball objects
  private final Point3f p1; // center of first endcap
  private final Point3f p2; // center of second endcap

  //TODO: store Ball references instead of points?
  //then we could extract new color information automatically
  private final Ball ball1; //Figure for first endpoint
  private final Ball ball2; //Figure for second endpoint
  
  //TODO: don't draw when balls on either end are not drawn
  //TODO: turning off bonds should set drawable to false; different from above
  //TODO: 'a' can cycle between cylinders, wireframe, and bonds off entirely
  private final Point3i p1i, p2i; // same in pixels, not shared
  private int bondType = CYLINDER;
protected float _d;
  
  public Bond(G3DSys g, Ball b0, Ball b1) {
    super(g);
    this.p1 = b0.getPoint();
    this.p2 = b1.getPoint();
    ball1 = b0;
    ball2 = b1;
    p1i = new Point3i();
    p2i = new Point3i();
  }
  
  public void draw() {
    //if atoms aren't drawn (because of atom filter), bonds aren't drawn
    if(!ball1.drawme || !ball2.drawme) return;
    //if points too far apart, assume molecule has been wrapped, ignore
    if(p1.distance(p2) > _gsys.getAngstromWidth()/2.0f) return;
    _gsys.screenSpace(p1, p1i);
    _gsys.screenSpace(p2, p2i);
    switch(bondType) {
    case CYLINDER:
      _gsys.getG3D().fillCylinder(ball1.getColor(),ball2.getColor(),
          Graphics3D.ENDCAPS_FLAT,
          _gsys.perspective(((int)(ball1.getPoint().z + ball2.getPoint().z))/2, getD()),
          p1i.x,p1i.y,p1i.z,p2i.x,p2i.y,p2i.z);
      break;
    case WIREFRAME:
      _gsys.getG3D().drawDashedLine(ball1.getColor(), 0, 0,
          p1i.x, p1i.y, p1i.z, p2i.x, p2i.y, p2i.z);
      //_gsys.getG3D().drawLine(color1,color2,p1i.x,p1i.y,p1i.z,p2i.x,p2i.y,p2i.z);
      break;
    }

  }

  public float getD() {
    // no good way to pick diameter yet
    // also, what does this mean when wireframe is on?
    return (ball1._d + ball2._d)*0.2f;
  }
  
  public void setBondType(int i) { bondType = i; }
  public int getBondType() { return bondType; }

  public Point3f getEndpoint1() { return p1; }
  public Point3f getEndpoint2() { return p2; }
  
  public short getColor1() { return ball1.getColor(); }
  public short getColor2() { return ball2.getColor(); }

/** set the molspace diameter of the figure, when applicable */
public void setD(float d) { _d = d; }
  
}
