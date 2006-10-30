package g3dsys.images;

import g3dsys.control.G3DSys;

public class TimedBall extends Ball {

	public TimedBall(G3DSys g, short c, float x, float y, float z, float diameter) {
		super(g, c, x, y, z, diameter);
	}
	
	public void draw() {
		/* An interesting result: ball rendering takes no less time when the window
		 * is clipped, but the drawing as a whole takes much shorter. Why?
		 */

		long start = System.currentTimeMillis();
		super.draw();
		//System.out.println("Render time for ball: "+(System.currentTimeMillis()-start));
	}
	
}
