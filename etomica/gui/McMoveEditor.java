package etomica.gui;
import etomica.*;
import java.beans.*;
import java.awt.*;

/**
 * Editor for an array of MCMoves (Monte Carlo moves).  
 * Provides a custom editor panel that permits selection of any instantiable McMove.
 */
public class McMoveEditor extends PropertyEditorSupport implements java.io.Serializable {
    
    private MCMove[] moves;
    IntegratorMC integrator;
    JavaWriter javaWriter;
    
    public McMoveEditor(IntegratorMC integrator) {
        this(integrator, null);
    }
    public McMoveEditor(IntegratorMC integrator, JavaWriter javaWriter) {
        this.integrator = integrator;
        this.javaWriter = javaWriter;
    }
    
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
    public void setValue(Object obj) {
        if(obj == null) return;
        MCMove[] newMoves = (MCMove[])obj;
        if(moves == null) { //no moves yet added
            moves = newMoves;
        }
        else { //adding to existing moves
            MCMove[] allMoves = new MCMove[moves.length + newMoves.length];
            for(int i=0; i<moves.length; i++) {
                allMoves[i] = moves[i];
            }
            for(int i=0; i<newMoves.length; i++) {
                allMoves[i+moves.length] = newMoves[i];
                if(javaWriter != null) javaWriter.add(newMoves[i]);
            }
            moves = allMoves;
        }
        firePropertyChange();
    }
    
    /**
     * Returns null to indicate no writing of initialization string.
     * Special method in JavaWriter for MCMove is used instead.
     */
    public String getJavaInitializationString() {return null;}
}