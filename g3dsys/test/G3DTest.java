package g3dsys.test;

import java.awt.Frame;
import java.awt.Panel;

import org.jmol.g3d.Graphics3D;

import g3dsys.control.G3DSys;
import junit.framework.TestCase;

public class G3DTest extends TestCase {
	
	private Frame frame;
	private Panel panel;
	private G3DSys gsys;
	
	public void setUp() {
		frame = new Frame();
		panel = new Panel();
		frame.add(panel);
		frame.setSize(400,400);
		panel.setSize(400,400);
		gsys = new G3DSys(panel);
	}
	
	/**
	 * Model should modify the Angstrom size as Figures are added and removed.
	 * Tests whether adding many randomly placed Figures and then removing
	 * them correctly leaves the size unchanged.
	 */
	public void testfigureRemoveResize() {
		gsys.addFig(G3DSys.BALL, Graphics3D.RED, 0,0,0, 100);
		float minx = gsys.getMinX(); float maxx = gsys.getMaxX();
		float miny = gsys.getMinY(); float maxy = gsys.getMaxY();
		float minz = gsys.getMinZ(); float maxz = gsys.getMaxZ();
		
		java.util.Random r = new java.util.Random();
		
		long begin = System.currentTimeMillis();
		java.util.Collection c = new java.util.ArrayList();
		for(int i=0; i<100000; i++) {
			float x = r.nextFloat()*r.nextInt();
			float y = r.nextFloat()*r.nextInt();
			float z = r.nextFloat()*r.nextInt();
			c.add(new Long( gsys.addFigNoRescale(G3DSys.BALL,Graphics3D.RED,x,y,z,100)));
		}
		gsys.recalcPPA();
		System.out.println("added figures in "+(System.currentTimeMillis()-begin)+" milliseconds");
		
		begin = System.currentTimeMillis();
		for(java.util.Iterator i = c.iterator(); i.hasNext();) {
			Long l = (Long)i.next();
			gsys.removeFigNoRescale(l);
		}
		gsys.removeFig(-1);
		System.out.println("removed figures in "+(System.currentTimeMillis()-begin)+" milliseconds");
		
		assertTrue("model not minx-resized properly; is "
				+gsys.getMinX()+", not "+minx, gsys.getMinX() == minx);
		assertTrue("model not maxx-resized properly; is "
				+gsys.getMaxX()+", not "+maxx, gsys.getMaxX() == maxx);
		assertTrue("model not miny-resized properly; is "
				+gsys.getMinY()+", not "+miny, gsys.getMinY() == miny);
		assertTrue("model not maxy-resized properly; is "
				+gsys.getMaxY()+", not "+maxy, gsys.getMaxY() == maxy);
		assertTrue("model not minz-resized properly; is "
				+gsys.getMinZ()+", not "+minz, gsys.getMinZ() == minz);
		assertTrue("model not maxz-resized properly; is "
				+gsys.getMaxZ()+", not "+maxz, gsys.getMaxZ() == maxz);
	}

}
