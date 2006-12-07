package etomica.graphics;
import java.awt.Color;
import java.util.HashMap;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

/**
 * Colors the atom according to the color given by its type field.
 *
 * @author David Kofke
 */
public final class ColorSchemeByType extends ColorScheme {
    
    public ColorSchemeByType() {
        colorMap = new HashMap();
    }
  
    public void setColor(AtomType type, Color c) {
        colorMap.put(type,c);
    }
    
    public Color getAtomColor(AtomLeaf a) {
        return getColor(a.getType());
    }
    
    public Color getColor(AtomType type) {
        Color color = (Color)colorMap.get(type);
        if (color == null) {
            color = defaultColor;
        }
        return color;
    }
    
    private final HashMap colorMap;
}
