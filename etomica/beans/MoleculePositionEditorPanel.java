package simulate;
import simulate.*;
import java.awt.*;
import java.beans.PropertyEditorSupport;
import javax.swing.JPanel;

public class MoleculePositionEditorPanel extends JPanel {

    DisplayPhase displayPhase;
    Phase phase;
    PropertyEditorSupport editor;
    Space.Coordinate[] coordinates;
    
    public MoleculePositionEditorPanel(PropertyEditorSupport ed) {
        editor = ed;
        setArray((Space.Coordinate[])ed.getValue());
        
        displayPhase = new DisplayPhase();
        this.add(displayPhase.graphic(null));
        
        displayPhase.setPhase(phase);
    }
    
    public void setArray(Space.Coordinate[] c) {
        if(c == null) coordinates = new Space.Coordinate[0];
        else {
            coordinates = c;
            if(c.length > 0) phase = c[0].parent().parentPhase();
        }
    }
    
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}