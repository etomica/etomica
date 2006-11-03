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

	public Display(G3DSys g) {
		gsys = g;
	}
	
	public void paint(Graphics g) {
		Graphics3D g3d = gsys.getG3D();
//		g3d.beginRendering(0, 0, gsys.getPixelWidth(),
//				gsys.getPixelHeight(), m, false);
		
		g3d.beginRendering(0,0,0,0, m, false); // ?!?!
		gsys.draw(); //dispatch to FigureManager
		g3d.endRendering();
		g.drawImage(g3d.getScreenImage(), 0, 0, null);
		g3d.releaseScreenImage();

	}
	
	public void update(Graphics g) {
		paint(g);
	}

}
