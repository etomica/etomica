
// Support for drawing a property value in a Canvas.

package etomica.gui;

import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.JInternalFrame;
import javax.swing.JPanel;

class PropertyCanvas extends JPanel implements MouseListener {

    PropertyCanvas(JInternalFrame frame, PropertyEditor pe) {
	    this.frame = frame;
	    editor = pe;
	    addMouseListener(this);
    }

    public void paint(Graphics g) {
	    Rectangle box = new Rectangle(2, 2, getSize().width - 4, getSize().height - 4);
	    editor.paintValue(g, box);
    }

    private static boolean ignoreClick = false;

    public void mouseClicked(MouseEvent evt) {
	    if (! ignoreClick) {
	        try {
		        ignoreClick = true;
		        int x = frame.getLocation().x - 30;
		        int y = frame.getLocation().y + 50;
		        new PropertyDialog(etomica.gui.Etomica.DesktopFrame.etomicaFrame, editor, x, y);
	        } 
	        finally {
    		    ignoreClick = false;
	        }
	    }
    }

    public void mousePressed(MouseEvent evt) {}

    public void mouseReleased(MouseEvent evt) {}

    public void mouseEntered(MouseEvent evt) {}

    public void mouseExited(MouseEvent evt) {}

    public String toString(){ return editor.getValue().getClass().getName(); }

    private JInternalFrame frame;
    public transient PropertyEditor editor;
}
