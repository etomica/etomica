package etomica;
import etomica.*;
import java.beans.*;
import java.awt.*;

/**
 * Editor for a property of class Meter[].  
 * Provides a custom editor panel that permits selection of any or all meters
 * from the list of those already instantiated.
 */
public class MeterArrayEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    private Meter[] meters;
    
    public String getAsText() {return null;}
    
    public boolean isPaintable() {return true;}
    
    public void paintValue(Graphics g, Rectangle box) {
        FontMetrics fm = g.getFontMetrics();
        String s = "Click to edit";
        int w = fm.stringWidth(s);
        int x = box.x;
        if(w < box.width) x += (box.width - w) / 2;
        int y = box.y + (box.height - fm.getHeight()) / 2 + fm.getAscent();
        g.drawString(s, x, y);
    }
    
    public boolean supportsCustomEditor() {return true;}
    
    public Component getCustomEditor() {
        return new MeterArrayEditorPanel(this);
    }
    
    public Object getValue() {return meters;}
    public void setValue(Object obj) {meters = (Meter[])obj;}
    
}