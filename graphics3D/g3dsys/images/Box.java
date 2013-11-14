package g3dsys.images;

import g3dsys.control.G3DSys;

import org.jmol.g3d.Graphics3D;
import org.jmol.util.Point3f;
import org.jmol.util.Point3i;


/**
 * Class for drawing a bounding box around the figures in molecule space
 */

public class Box extends Figure {
	
	// eight points to define the box
	private Point3f LUT,LDT,RUT,RDT; //leftUpTop, leftDownTop, etc.
	private Point3f LUB,LDB,RUB,RDB; //leftUpBottom, leftDownBottom, etc.
    private final Point3i t1, t2;
    protected short _c;

	public Box(G3DSys g, short c) {
		super(g);
		_c = c;
		resize();
		
        
        t1 = Point3i.new3(0,0,0);
        t2 = Point3i.new3(0,0,0);
		
        Point3f p = Point3f.new3(
                (_gsys.getMinX() + _gsys.getMaxX())/2,
                (_gsys.getMinY() + _gsys.getMaxY())/2,
                (_gsys.getMinZ() + _gsys.getMaxZ())/2);
		g.setCenterOfRotation(p); //rotate about center of box
	}

	private void resize() {
		float minX = _gsys.getMinX(); float maxX = _gsys.getMaxX();
		float minY = _gsys.getMinY(); float maxY = _gsys.getMaxY();
		float minZ = _gsys.getMinZ(); float maxZ = _gsys.getMaxZ();

		LUT = Point3f.new3(minX,minY,minZ); LUB = Point3f.new3(minX,minY,maxZ);
		LDT = Point3f.new3(minX,maxY,minZ); LDB = Point3f.new3(minX,maxY,maxZ);
		RUT = Point3f.new3(maxX,minY,minZ); RUB = Point3f.new3(maxX,minY,maxZ);
		RDT = Point3f.new3(maxX,maxY,minZ); RDB = Point3f.new3(maxX,maxY,maxZ);
	}

	public void draw() {
		Graphics3D g3d = _gsys.getG3D();
		
		if (!g3d.setColix(_c)) return;
        resize(); //must resize in case the model has changed

        //four top to bottom lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(LUB, t2);
		g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(LDT, t1); _gsys.screenSpace(LDB, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(RUT, t1); _gsys.screenSpace(RUB, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(RDT, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLineAB(t1, t2);
		
		//four left to right lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(RUT, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(LDT, t1); _gsys.screenSpace(RDT, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(LUB, t1); _gsys.screenSpace(RUB, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(LDB, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLineAB(t1, t2);
		
		//four up to down lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(LDT, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(RUT, t1); _gsys.screenSpace(RDT, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(LUB, t1); _gsys.screenSpace(LDB, t2);
        g3d.drawLineAB(t1, t2);
        _gsys.screenSpace(RUB, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLineAB(t1, t2);
	}
	
	public float getD() {
		return 0;
	}

    /** @return the color of the figure */
    public short getColor() { return _c; }

    /** set the color of the figure */
    public void setColor(short c) { _c = c; }

}
