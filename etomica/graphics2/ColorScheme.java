package etomica.graphics2;

import etomica.atom.IAtom;

/**
 * Defines an interface to retrieve the Atom color from an Atom object
 * @author Henrique Bucher
 */
 
public interface ColorScheme {
	public int getNumColors();
    public int atomColor(IAtom a);
    public Color getColor( int index );
}//end of ColorScheme
