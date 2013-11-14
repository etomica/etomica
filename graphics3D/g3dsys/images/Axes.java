package g3dsys.images;

import g3dsys.control.G3DSys;

import org.jmol.util.Point3f;
import org.jmol.util.Point3i;

public class Axes extends Figure {

	private final G3DSys gsys;
    private final Point3i t1, t2;
    protected short _c;
	
	public Axes(G3DSys g, short c) {
		super(g);
        _c = c;
		gsys = g;
        t1 = new Point3i();
        t2 = new Point3i();
	}

	public void draw() {
        if (!_gsys.getG3D().setColix(_c)) return;
		Point3f cor = gsys.getCenterOfRotation();
		Point3f xNeg = Point3f.new3(gsys.getMinX(),cor.y,cor.z);
		Point3f xPos = Point3f.new3(gsys.getMaxX(),cor.y,cor.z);
		Point3f yNeg = Point3f.new3(cor.x,gsys.getMinY(),cor.z);
		Point3f yPos = Point3f.new3(cor.x,gsys.getMaxY(),cor.z);
		Point3f zNeg = Point3f.new3(cor.x,cor.y,gsys.getMinZ());
		Point3f zPos = Point3f.new3(cor.x,cor.y,gsys.getMaxZ());
        _gsys.screenSpace(xNeg, t1); _gsys.screenSpace(xPos, t2);
        _gsys.getG3D().drawLineAB(t1, t2);
        _gsys.screenSpace(yNeg, t1); _gsys.screenSpace(yPos, t2);
        _gsys.getG3D().drawLineAB(t1, t2);
        _gsys.screenSpace(zNeg, t1); _gsys.screenSpace(zPos, t2);
        _gsys.getG3D().drawLineAB(t1, t2);
	}

	public float getD() { return 0; }

    /** @return the color of the figure */
    public short getColor() { return _c; }

    /** set the color of the figure */
    public void setColor(short c) { _c = c; }

}
