package etomica.graphics2;

import etomica.api.IVector;

/** Interface used by SceneManager to communicate with the graphics package. 
 *  Should be implemented by all classes willing to draw the contents of etomica. */

public interface Renderable
{

	public interface Shape 
	{
		public void setPosition( IVector pos );
		public void setScale( IVector scale );
		public void setRotation( IVector fromvec, IVector tovec );
		public void setRotation( float angle, IVector axis );
		public void setColor( int color_index );
		public void setColorScheme( ColorScheme scheme );
        public void dispose();
	};
    
    public void setColorScheme(ColorScheme scheme);
    
	public interface Sphere extends Shape
	{
	};
	public interface Polyline extends Shape
	{
		public void appendLine( IVector pointA, IVector pointB );
	};
	public Sphere createSphere();
	public Polyline createPoly();
	
};
