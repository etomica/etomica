package etomica;
import etomica.*;
import java.beans.*;
import java.awt.*;

/**
 * Editor for an array of MCMoves (Monte Caro moves).  
 * Provides a custom editor panel that permits selection of any instantiable McMove.
 */
public class McMoveEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    private MCMove[] moves;
    
    public String getAsText() {return null;}
    
    public boolean isPaintable() {return true;}
    
    public void paintValue(Graphics g, Rectangle box) {
        FontMetrics fm = g.getFontMetrics();
        String s = "Double click to edit";
        int w = fm.stringWidth(s);
        int x = box.x;
        if(w < box.width) x += (box.width - w) / 2;
        int y = box.y + (box.height - fm.getHeight()) / 2 + fm.getAscent();
        g.drawString(s, x, y);
    }
    
    public boolean supportsCustomEditor() {return true;}
    
    public Component getCustomEditor() {
        return new McMoveEditorPanel(this);
    }
    
    public Object getValue() {return moves;}
    public void setValue(Object obj) {moves = (MCMove[])obj;}
    
 /*   public String getJavaInitializationString() {
        System.out.println("inside ModulatorEditor.getJavaInitializationString");
        if(getValue() == null) return "";
        else {
            String[] value = (String[])getValue();
            return "new Modulator("+value[0]+",\""+value[1]+"\")";
        }
    }*/
}