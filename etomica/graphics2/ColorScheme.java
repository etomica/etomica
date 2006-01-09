package etomica.graphics2;

import etomica.atom.AtomLeaf;

/**
 * Defines an interface to retrieve the Atom color from an Atom object
 * @author Henrique Bucher
 */
 
public interface ColorScheme {
	public int getNumColors();
    public int atomColor(AtomLeaf a);
    public Color getColor( int index );
}//end of ColorScheme
