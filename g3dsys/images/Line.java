package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;


public class Line extends Figure {
	
	protected Point3f start, end;
	
	public Line(G3DSys m, short c, Point3f s, Point3f e) {
		super(m,c);
		start = s;
		end = e;
	}

	public void draw() {
		Point3i s = _gsys.screenSpace(start);
		Point3i e = _gsys.screenSpace(end);
		_gsys.getG3D().drawDashedLine(_c, 0, 0, s, e);
	}

	public float getD() { return 0; }
	public void setEnd(Point3f p) { end = p; }

}
