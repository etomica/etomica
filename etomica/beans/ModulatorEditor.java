package etomica;
import etomica.*;
import java.beans.*;
import java.awt.*;

/**
 * Editor for a property of class Modulator.  
 * Provides a custom editor panel that permits selection of any object in
 * the form designer and then any property of the selected object.  Selection
 * results in a JavaInitializationString the prescribes instantiation of 
 * a new modulator for the selected object/property
 */
public class ModulatorEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    private ModulatorAbstract modulator;
    
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
        return new ModulatorEditorPanel(this);
    }
    
    public Object getValue() {return modulator;}
    public void setValue(Object obj) {modulator = (ModulatorAbstract)obj;}
    
    public String getJavaInitializationString() {
        System.out.println("inside ModulatorEditor.getJavaInitializationString");
        if(getValue() == null) return "";
        else {
            String[] value = (String[])getValue();
            return "new Modulator("+value[0]+",\""+value[1]+"\")";
        }
    }
}