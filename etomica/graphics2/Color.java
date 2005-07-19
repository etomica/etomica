/*
 * Created on Jul 19, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.graphics2;

/**
 * @author Henrique
 * A class to represent the color in a device-independent way (not dependent on AWT or SWT)
 */
public class Color {
	public Color( float cr, float cg, float cb ) { r=cr; g=cg; b=cb; }
	public float r,g,b;
	public static final Color RED = new Color(1.0f,0,0);
	public static final Color BLUE = new Color(0,0,1.0f);
	public static final Color GREEN = new Color(0,1.0f,0);
	public static final Color WHITE = new Color(1,1,1);
	public static final Color BLACK = new Color(0,0,0);
	public static final Color GRAY75 = new Color(0.75f,0.75f,0.75f);
	public static final Color GRAY50 = new Color(0.5f,0.5f,0.5f);
	public static final Color GRAY25 = new Color(0.25f,0.25f,0.25f);
}
