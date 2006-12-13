package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

public class Axes extends Figure {

	private final G3DSys gsys;
    private final Point3i t1, t2;
	
	public Axes(G3DSys g, short c) {
		super(g, c);
		gsys = g;
        t1 = new Point3i();
        t2 = new Point3i();
	}

	public void draw() {
		Point3f cor = gsys.getCenterOfRotation();
		Point3f xNeg = new Point3f(gsys.getMinX(),cor.y,cor.z);
		Point3f xPos = new Point3f(gsys.getMaxX(),cor.y,cor.z);
		Point3f yNeg = new Point3f(cor.x,gsys.getMinY(),cor.z);
		Point3f yPos = new Point3f(cor.x,gsys.getMaxY(),cor.z);
		Point3f zNeg = new Point3f(cor.x,cor.y,gsys.getMinZ());
		Point3f zPos = new Point3f(cor.x,cor.y,gsys.getMaxZ());
        _gsys.screenSpace(xNeg, t1); _gsys.screenSpace(xPos, t2);
        _gsys.getG3D().drawLine(_c, t1, t2);
        _gsys.screenSpace(yNeg, t1); _gsys.screenSpace(yPos, t2);
        _gsys.getG3D().drawLine(_c, t1, t2);
        _gsys.screenSpace(zNeg, t1); _gsys.screenSpace(zPos, t2);
        _gsys.getG3D().drawLine(_c, t1, t2);
	}

	public float getD() { return 0; }

}
