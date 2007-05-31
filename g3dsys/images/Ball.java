package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

/**
 * Class for drawing spheres (atoms)
 */

public class Ball extends Figure {
	
	private Point3f p; // xyz position in molecule space
    private Point3i s;
    protected short _c;
    protected Point3f _p;
    protected float _d;
    
	public Ball(g3dsys.control.G3DSys g, short c,
			float x, float y, float z, float diameter) {
		this(g, c, new Point3f(x,y,z), diameter);
	}
	public Ball(g3dsys.control.G3DSys g, short c, Point3f point, float diameter) {
		super(g);
        _c = c;
		_d = diameter;
		p = point;
		if(p != null) _p = p;
        s = new Point3i();
	}
	
	public void draw() {
        if(!drawme) return;
		//if overhead is too much, give figures actual references later
        _gsys.screenSpace(p, s);
		int diam = _gsys.perspective(s.z, _d);
		_gsys.getG3D().fillSphereCentered(_c, diam, s);
	}
	
	public float getD() { return _d; }
    /** @return the color of the figure */
    public short getColor() { return _c; }
    /** set the color of the figure */
    public void setColor(short c) { _c = c; }
    /** @return the molspace position of the figure */
    public Point3f getPoint() { return _p; }
    /** @return the x molspace position of the figure */
    public float getX() { return _p.x; }
    /** set the x molspace position of the figure */
    public void setX(float x) { _p.x = x; }
    /** @return the y molspace position of the figure */
    public float getY() { return _p.y; }
    /** set the y molspace position of the figure */
    public void setY(float y) { _p.y = y; }
    /** @return the z molspace position of the figure */
    public float getZ() { return _p.z; }
    /** set the z molspace position of the figure */
    public void setZ(float z) { _p.z = z; }
    /** set the molspace diameter of the figure, when applicable */
    public void setD(float d) { _d = d; }

}
