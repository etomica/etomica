package g3dsys.images;

import org.jmol.util.Point3f;
import org.jmol.util.Point3i;

/**
 * Class for drawing ellipses
 */

public class Ellipse extends Figure {
	
	private Point3f p; // xyz position in molecule space
    private Point3i s;
    protected short _c;
    protected Point3f _p;
    protected float _d;
    protected float _a;
    protected final Ellipse3D myEllipse3D;
    
	public Ellipse(g3dsys.control.G3DSys g, short c,
			float x, float y, float z, float diameter) {
		this(g, c, Point3f.new3(x,y,z), diameter, 1);
	}
	public Ellipse(g3dsys.control.G3DSys g, short c, Point3f point, float diameter, float aspect) {
		super(g);
        _c = c;
		_d = diameter;
        _a = aspect;
		p = point;
		if(p != null) _p = p;
        s = new Point3i();
        myEllipse3D = new Ellipse3D(g.getG3D());
	}
	
	public void draw() {
        if(!drawme) return;
        if (!_gsys.getG3D().setColix(_c)) return;
		//if overhead is too much, give figures actual references later
        _gsys.screenSpace(p, s);
		int diam = _gsys.perspective(s.z, _d);
		myEllipse3D.plotCircleCenteredClipped(s.x, s.y, s.z, diam, _a);
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
    /** set the aspect ratio (y/x) of the figure */
    public float getAspect() { return _a; }
    /** @return the aspect ratio (y/x) of the figure */
    public void setAspect(float a) { _a = a; }

}
