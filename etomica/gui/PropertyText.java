
// Support for a PropertyEditor that uses text.

package etomica.gui;

import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.JTextField;

public class PropertyText extends JTextField implements KeyListener, FocusListener, PropertyChangeListener {

    public PropertyText(PropertyEditor pe) {
	    super(pe.getAsText());
	    editor = pe;
	    addKeyListener(this);
	    addFocusListener(this);
	    editor.addPropertyChangeListener(this);
    	setBorder(PropertySheet.EMPTY_BORDER);
    }

    public void repaint() {}

    protected void updateEditor() {
	    try {
	        editor.setAsText(getText());
	    } 
	    catch (IllegalArgumentException ex) {
	        // Quietly ignore.
	    }
    }
    
    /**
     * Listen to update display if editor changes value in some other way.
     * For example, display of dimensioned property values can be changed with a change of the units.
     */
    public void propertyChange(PropertyChangeEvent evt) {
        setText(editor.getAsText());
    }
    
    //----------------------------------------------------------------------
    // Focus listener methods.

    public void focusGained(FocusEvent e) {}

    public void focusLost(FocusEvent e) {
    	updateEditor();
    }
    
    //----------------------------------------------------------------------
    // Keyboard listener methods.

    public void keyReleased(KeyEvent e) {
 	    if (e.getKeyCode() == KeyEvent.VK_ENTER) {
	        updateEditor();
	    }
    }

    public void keyPressed(KeyEvent e) {}

    public void keyTyped(KeyEvent e) {}

    //----------------------------------------------------------------------
    private transient PropertyEditor editor;
}
