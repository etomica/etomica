package etomica.modules.crystalviewer;
import java.awt.Color;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.ISimulation;
import etomica.graphics.ColorScheme;
import etomica.lattice.LatticePlane;

/**
 * Colors atoms according to whether or not they are in a plane,
 * as determined by the inPlane method of a Plane instance.
 *
 * @author David Kofke
 * @see etomica.lattice.LatticePlane;
 */
public class ColorSchemePlane extends ColorScheme {
    
    private Color colorIn, colorOut;
    private LatticePlane plane;
    
    public ColorSchemePlane(ISimulation sim, LatticePlane plane) {this(sim, plane, Color.red, Color.green);}
    public ColorSchemePlane(ISimulation sim, LatticePlane plane, Color colorIn, Color colorOut) {
        super(sim, colorOut);
        this.plane = plane;
        this.colorIn = colorIn;
        this.colorOut = colorOut;
    }
    public Color getAtomColor(IAtomLeaf a) {
        return plane.inPlane(((IAtomPositioned)a).getPosition()) ? colorIn : colorOut;
    }
    
    /**
     * Color for atoms that are in the plane.
     */
    public void setColorIn(Color colorIn) {this.colorIn = colorIn;}
    public Color getColorIn() {return colorIn;}
    
    /**
     * Color for atoms that are out of the plane.
     */
    public void setColorOut(Color colorOut) {this.colorOut = colorOut;}
    public Color getColorOut() {return colorOut;}
}//end of Simple
