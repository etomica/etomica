package etomica.graphics2;

import etomica.Atom;
import etomica.AtomType;
import java.util.HashMap;
import java.util.LinkedList;

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
    public final int atomColor(Atom a) 
    {
    	Integer index = (Integer) colormap.get( a.type );
    	return index.intValue();
    }

    public void setColor(AtomType type, Color c) 
    {
    	int index = colorindex.size();
    	colormap.put( type, new Integer(index) );
    	colorindex.add( c );
    }


	protected HashMap colormap = new HashMap();
	protected LinkedList colorindex = new LinkedList();
	
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
		return (Color)colorindex.get( index );
	}
   
}
