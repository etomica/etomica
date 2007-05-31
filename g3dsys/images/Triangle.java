package g3dsys.images;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

/**
 * Class for drawing triangles (possibly making up a polygon)
 */

public class Triangle extends Figure {
	
	private Point3f p1, p2, p3;
    private Point3i s1, s2, s3;
    protected short _c;
    
	public Triangle(g3dsys.control.G3DSys g, short c,
			Point3f vertex1, Point3f vertex2, Point3f vertex3) {
		super(g);
        _c = c;
        p1 = new Point3f(vertex1);
        p2 = new Point3f(vertex2);
        p3 = new Point3f(vertex3);
        s1 = new Point3i();
        s2 = new Point3i();
        s3 = new Point3i();
	}
	
	public void draw() {
        if(!drawme) return;
		//if overhead is too much, give figures actual references later
        _gsys.screenSpace(p1, s1);
        _gsys.screenSpace(p2, s2);
        _gsys.screenSpace(p3, s3);
		_gsys.getG3D().fillTriangle(_c, s1, s2, s3);
	}
	
    /** @return the color of the figure */
    public short getColor() { return _c; }
    /** set the color of the figure */
    public void setColor(short c) { _c = c; }

    public Point3f getVertex1() {
        return p1;
    }
    public Point3f getVertex2() {
        return p2;
    }
    public Point3f getVertex3() {
        return p3;
    }
}
