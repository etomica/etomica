package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

/**
 * Class for drawing spheres (atoms)
 */

public class Ball extends Figure {
	
	private Point3f p; // xyz position in molecule space
	
	public Ball(g3dsys.control.G3DSys g, short c,
			float x, float y, float z, float diameter) {
		this(g, c, new Point3f(x,y,z), diameter);
	}
	public Ball(g3dsys.control.G3DSys g, short c, Point3f point, float diameter) {
		super(g,c);
		_d = diameter;
		p = point;
		if(p != null) _p = p;
	}
	
	public void draw() {
        if(!drawme) return;
		//if overhead is too much, give figures actual references later
		Point3i s = _gsys.screenSpace(p);
		int diam = _gsys.perspective(s.z, _d);
		_gsys.getG3D().fillSphereCentered(_c, diam, s);
	}
	
	public float getD() { return _d; }

}
