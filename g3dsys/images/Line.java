package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;


public class Line extends Figure {
	
	protected Point3f start, end;
    protected Point3i s, e;
	
    public Line(G3DSys m, short c, Point3f start, Point3f end) {
        super(m,c);
        this.start = start;
        this.end = end;
        s = new Point3i();
        e = new Point3i();
    }

    public void draw() {
        _gsys.screenSpace(start, s);
        _gsys.screenSpace(end, e);
        _gsys.getG3D().drawDashedLine(_c, 0, 0, s, e);
    }

	public float getD() { return 0; }
	public void setEnd(Point3f p) { end = p; }

}
