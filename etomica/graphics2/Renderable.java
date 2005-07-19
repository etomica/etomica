package etomica.graphics2;

/** Interface used by SceneManager to communicate with the graphics package. 
 *  Should be implemented by all classes willing to draw the contents of etomica. */

public interface Renderable
{
	// SHAPE TYPES
	public final int ELLIPSE = 1;
	
	// PROPERTY IDS
	public final int ELLIPSE_RADIUS = 1001; // float=>sphere, float[3]=>ellipse
	
	public int  createObject( int type );
	public void setObjectPosition( int index, float x, float y, float z );
	public void setObjectColor( int index, int color_index );
	public void setObjectProperty( int index, int prop_type, float value );
	public void setObjectProperty( int index, int prop_type, float[] value );
	public void setObjectProperty( int index, int prop_type, Object value );
	
	public void setColorScheme( ColorScheme scheme );
	
	public float[] getCameraPosition();
};
