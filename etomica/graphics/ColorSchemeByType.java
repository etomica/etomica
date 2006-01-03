package etomica.graphics;
import java.awt.Color;
import java.util.HashMap;

import etomica.atom.Atom;
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
    
    public Color getAtomColor(Atom a) {
        return getColor(a.type);
    }
    
    public Color getColor(AtomType type) {
        Color color = (Color)colorMap.get(type);
        if(color == null) color = defaultColor;
        return color;
    }
    
   /**
     * @return Returns the defaultColor.
     */
    public Color getDefaultColor() {
        return defaultColor;
    }
    /**
     * @param defaultColor The defaultColor to set.
     */
    public void setDefaultColor(Color defaultColor) {
        this.defaultColor = defaultColor;
    }
    
    private Color defaultColor = Color.RED;
    private final HashMap colorMap;
}
