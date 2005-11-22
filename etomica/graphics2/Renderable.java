package etomica.graphics2;

import etomica.space.Vector;

/** Interface used by SceneManager to communicate with the graphics package. 
 *  Should be implemented by all classes willing to draw the contents of etomica. */

public interface Renderable
{

	public interface Shape 
	{
		public void setPosition( Vector pos );
		public void setScale( Vector scale );
		public void setRotation( Vector fromvec, Vector tovec );
		public void setRotation( float angle, Vector axis );
		public void setColor( int color_index );
		public void setColorScheme( ColorScheme scheme );
        public void dispose();
	};
	public interface Sphere extends Shape
	{
	};
	public interface Polyline extends Shape
	{
		public void appendLine( Vector pointA, Vector pointB );
	};
	public Sphere createSphere();
	public Polyline createPoly();
	
};
