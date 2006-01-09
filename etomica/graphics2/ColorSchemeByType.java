package etomica.graphics2;

import java.util.ArrayList;
import java.util.HashMap;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

/**
 * Colors the atom according to the color given by its type field.  
 *
 * @author David Kofke
 * @author Henrique Bucher
 */

public final class ColorSchemeByType 
implements ColorScheme, java.io.Serializable 
{
    
	public ColorSchemeByType() 	{}

 /**
  * Initialize atom color to the color of its type
  */
    public final int atomColor(AtomLeaf a) 
    {
    	Integer index = (Integer) colormap.get( a.type );
    	if ( index==null )
    	{
    		// Assign the next position to this unknown
    		Color nextcolor = DEFAULTCOLORLIST[position_in_list++];
    		if ( position_in_list==DEFAULTCOLORLIST.length )
    			position_in_list = 0;
    		colormap.put( a.type, new Integer( colorindex.size() ) );
    		colorindex.add( nextcolor );
    		return atomColor( a );
    		
    	}
    	return index.intValue();
    }

    public void setColor(AtomType type, Color c) 
    {
    	int index = colorindex.size();
    	colormap.put( type, new Integer(index) );
    	colorindex.add( c );
    }


	protected HashMap colormap = new HashMap();
	protected ArrayList colorindex = new ArrayList();
	
	/** Get the total number of colors for this colorscheme
	 * @see etomica.graphics2.ColorScheme#getNumColors()
	 */
	public int getNumColors() {
		return colorindex.size();
	}

	/* (non-Javadoc)
	 * @see etomica.graphics2.ColorScheme#getColor(int)
	 */
	public Color getColor(int index) {
		if ( index==-1 )
			return DEFAULTCOLOR;
		return (Color)colorindex.get( index );
	}
	
	static private final Color[] DEFAULTCOLORLIST = new Color[]{ Color.RED, Color.GREEN, Color.BLUE, Color.WHITE, Color.BLACK, Color.GRAY25, Color.GRAY50, Color.GRAY75 };
	int position_in_list =0;
	static private final Color DEFAULTCOLOR = new Color( 0.5f, 0.5f, 0.5f );
   
}
