package etomica;
import etomica.*;
import java.beans.*;
import java.awt.*;

public class MoleculePositionEditor extends PropertyEditorSupport {
    
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
        return new MoleculePositionEditorPanel(this);
    }
    
    public void setValue(Object value) {
        if(value==null) System.out.println("null passed to MoleculePositionEditor.setValue()");
        else System.out.println("inside MoleculePositionEditor.setValue "+value.toString());
        super.setValue(value);
        
    }
    
    public String getJavaInitializationString() {
        System.out.println("inside MoleculePositionEditor.getJavaInitializationString");
        return "true";
    }
}