package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;


public class Line extends Figure {
	
  protected Point3f start, end;
  protected Point3i s, e;
protected short _c;
	
  public Line(G3DSys g, short c, Point3f start, Point3f end) {
    super(g);
    _c = c;
    this.start = start;
    this.end = end;
    s = new Point3i();
    e = new Point3i();
  }

  public void draw() {
    if(!drawme) return;
    _gsys.screenSpace(start, s);
    _gsys.screenSpace(end, e);
    _gsys.getG3D().drawDashedLine(_c, 0, 0, s.x, s.y, s.z, e.x, e.y, e.z);
  }

  public float getD() { return 0; }
  public void setEnd(Point3f p) { end = p; }

  public void setStart(float x, float y, float z) {
    start.set(x, y, z);
  }
    
  public void setEnd(float x, float y, float z) {
    end.set(x, y, z);
  }
    
  public Point3f getStart() { return start; }
  public Point3f getEnd() { return end; }

/** @return the color of the figure */
public short getColor() { return _c; }

/** set the color of the figure */
public void setColor(short c) { _c = c; }
}
