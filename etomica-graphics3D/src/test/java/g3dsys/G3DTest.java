/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package g3dsys;

import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Figure;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class G3DTest {
	
	private Frame frame;
	private Panel panel;
	private G3DSys gsys;
	
	@BeforeEach
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
	@Test
	public void testfigureRemoveResize() {
        gsys.addFig(new Ball(gsys, G3DSys.getColix(Color.RED), 0,0,0, 100));
		float minx = gsys.getMinX(); float maxx = gsys.getMaxX();
		float miny = gsys.getMinY(); float maxy = gsys.getMaxY();
		float minz = gsys.getMinZ(); float maxz = gsys.getMaxZ();
		
		java.util.Random r = new java.util.Random();
		
		long begin = System.currentTimeMillis();
		List<Ball> c = new ArrayList<Ball>();
		for(int i=0; i<100000; i++) {
			float x = r.nextFloat()*r.nextInt();
			float y = r.nextFloat()*r.nextInt();
			float z = r.nextFloat()*r.nextInt();
            Ball ball = new Ball(gsys,G3DSys.getColix(Color.RED),x,y,z,100);
            gsys.addFig(ball);
            c.add(ball);
		}
		System.out.println("added figures in "+(System.currentTimeMillis()-begin)+" milliseconds");
		
		begin = System.currentTimeMillis();
		for(Iterator<Ball> i = c.iterator(); i.hasNext();) {
            Figure f = i.next();
            gsys.removeFig(f);
        }
		System.out.println("removed figures in "+(System.currentTimeMillis()-begin)+" milliseconds");
		
		Assertions.assertTrue(gsys.getMinX() == minx, "model not minx-resized properly; is "
				+ gsys.getMinX() + ", not " + minx);
		Assertions.assertTrue(gsys.getMaxX() == maxx, "model not maxx-resized properly; is "
				+ gsys.getMaxX() + ", not " + maxx);
		Assertions.assertTrue(gsys.getMinY() == miny, "model not miny-resized properly; is "
				+ gsys.getMinY() + ", not " + miny);
		Assertions.assertTrue(gsys.getMaxY() == maxy, "model not maxy-resized properly; is "
				+ gsys.getMaxY() + ", not " + maxy);
		Assertions.assertTrue(gsys.getMinZ() == minz, "model not minz-resized properly; is "
				+ gsys.getMinZ() + ", not " + minz);
		Assertions.assertTrue(gsys.getMaxZ() == maxz, "model not maxz-resized properly; is "
				+ gsys.getMaxZ() + ", not " + maxz);
	}

}
