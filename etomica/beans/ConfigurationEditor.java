package simulate;
import simulate.*;
import java.beans.*;
import java.awt.*;

/**
 * Editor for a property of class Configuration.  
 * Provides a custom editor panel that permits selection of configurations
 *
 * @author David Kofke
 */
public class ConfigurationEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    private Configuration configuration;
    
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
        return new ConfigurationEditorPanel(this);
    }
    
    public Object getValue() {return configuration;}
    public void setValue(Object obj) {configuration = (Configuration)obj;}
    
}