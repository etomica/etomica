package g3dsys.control;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Panel;

import javax.vecmath.Matrix3f;

import org.jmol.g3d.Graphics3D;

/**
 *	Represents the panel on which all graphical data is actually drawn.
 *	Intercepts mouse events, etc. 
 *
 */

class Display extends Panel {
	
	private G3DSys gsys;
	private Matrix3f m = new Matrix3f(1,0,0,0,1,0,0,0,1);

	//double buffer
	private Graphics buf;
	private Image offscr;
	private Dimension size;
	
	public Display(G3DSys g) {
		gsys = g;
		
		size = getSize();
		offscr = createImage(size.width, size.height);
		//buf = offscr.getGraphics();
	}
	
	public void paint(Graphics g) {

		//seems to be using same clip even when partially offscreen
		//System.out.println("using clip: "+g.getClipBounds());

		Graphics3D g3d = gsys.getG3D();
//		g3d.beginRendering(0, 0, gsys.getPixelWidth(),
//				gsys.getPixelHeight(), m, false);
		
//		buf.clearRect(0,0,size.width,size.width);
//    buf.setColor(java.awt.Color.red);
//    buf.drawString("Bad Double-buffered",10,10);
		
		g3d.beginRendering(0,0,0,0, m, false); // ?!?!
		long start = System.currentTimeMillis();

		gsys.draw(); //dispatch to FigureManager

		g3d.endRendering();
		System.out.println(g3d.getScreenImage().getWidth(null));
		g.drawImage(g3d.getScreenImage(), 0, 0, null);
		System.out.println("END "+(System.currentTimeMillis()-start));		

		g3d.releaseScreenImage();

	}
	
	public void update(Graphics g) {
		paint(g);
	}

}
