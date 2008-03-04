package etomica.graphics;
import java.awt.Color;
import java.util.HashMap;

import etomica.api.IAtom;
import etomica.api.IAtomType;

/**
 * Colors the atom according to the color given by its type field.
 *
 * @author David Kofke
 */
public final class ColorSchemeByType extends ColorScheme {
    
    public ColorSchemeByType() {
        colorMap = new HashMap();
    }
  
    public void setColor(IAtomType type, Color c) {
        colorMap.put(type,c);
    }
    
    public Color getAtomColor(IAtom a) {
        return getColor(a.getType());
    }
    
    public Color getColor(IAtomType type) {
        Color color = (Color)colorMap.get(type);
        if (color == null) {
            if (defaultColorsUsed < moreDefaultColors.length) {
                color = moreDefaultColors[defaultColorsUsed];
                defaultColorsUsed++;
                setColor(type, color);
            }
            else {
                color = defaultColor;
                setColor(type, color);
            }
        }
        return color;
    }
    
    private final HashMap colorMap;
    protected final Color[] moreDefaultColors = new Color[]{ColorScheme.DEFAULT_ATOM_COLOR, Color.BLUE, Color.GREEN, Color.YELLOW, Color.ORANGE};
    protected int defaultColorsUsed = 0;
}
