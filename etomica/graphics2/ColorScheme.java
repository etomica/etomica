package etomica.graphics2;

import etomica.api.IAtomLeaf;

/**
 * Defines an interface to retrieve the Atom color from an Atom object
 * @author Henrique Bucher
 */
 
public interface ColorScheme {
	public int getNumColors();
    public int atomColor(IAtomLeaf a);
    public Color getColor( int index );
}//end of ColorScheme
