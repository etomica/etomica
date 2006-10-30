package g3dsys.images;

import javax.vecmath.Point3f;

import g3dsys.control.G3DSys;

public class Axes extends Figure {

	private G3DSys gsys;
	
	public Axes(G3DSys g, short c) {
		super(g, c);
		gsys = g;
	}

	public void draw() {
		Point3f cor = gsys.getCenterOfRotation();
		Point3f xNeg = new Point3f(gsys.getMinX(),cor.y,cor.z);
		Point3f xPos = new Point3f(gsys.getMaxX(),cor.y,cor.z);
		Point3f yNeg = new Point3f(cor.x,gsys.getMinY(),cor.z);
		Point3f yPos = new Point3f(cor.x,gsys.getMaxY(),cor.z);
		Point3f zNeg = new Point3f(cor.x,cor.y,gsys.getMinZ());
		Point3f zPos = new Point3f(cor.x,cor.y,gsys.getMaxZ());
		gsys.getG3D().drawLine(_c, gsys.screenSpace(xNeg), gsys.screenSpace(xPos));
		gsys.getG3D().drawLine(_c, gsys.screenSpace(yNeg), gsys.screenSpace(yPos));
		gsys.getG3D().drawLine(_c, gsys.screenSpace(zNeg), gsys.screenSpace(zPos));
	}

	public float getD() { return 0; }

}
