package etomica;
import etomica.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;

// attempt at a generic editor panel for an array of simulation elements (e.g., Meter[], Phase[]).
// doesn't work because of problem casting Object[] to a specific array type


/**
 * Panel for use with a Meter[] property editor.  Permits selection of
 * set of meters from those that are presently instantiated.  Presents
 * an array of checkboxes for the selection.
 */
public class ElementArrayEditorPanel extends JPanel {

    PropertyEditorSupport editor;
    MyCheckBox[] checkBox; //array of checkboxes with associated objects
    
    public ElementArrayEditorPanel(PropertyEditorSupport ed, DefaultListModel elementList) {
        editor = ed;
        
        //discover objects presently instantiated
        Object[] objects = new Object[elementList.size()];
        for (int i = 0; i < elementList.size(); i++)
            objects[i] = elementList.elementAt(i);
        
        //construct panel of checkboxes
        Object[] currentObjects = (Object[])editor.getValue();
        JPanel objectPanel = new JPanel(new GridLayout(0,2));
        checkBox = new MyCheckBox[objects.length];
        for(int i=0; i<objects.length; i++) {
            checkBox[i] = new MyCheckBox(objects[i]);
            objectPanel.add(checkBox[i]);
            for(int j=0; j<currentObjects.length; j++) {  //put check in box if already in list of objects
                if(currentObjects[j] == objects[i]) checkBox[i].setSelected(true);
            }
        }
        JButton button = new JButton("OK");
        button.addActionListener(new DoneButtonListener());
                
        add(button);
        add(objectPanel);
    }
    
    /**
     * Simple extension of JCheckBox to permit it to be associated with an object
     */
    private class MyCheckBox extends JCheckBox {
        Object object;
        MyCheckBox(Object obj) {
            super(obj.toString());
            object = obj;
        }
    }
    
	
    /**
     * ActionListener class called when a property is selected.
     * Action is to set value in editor and fire property change event.
     */
	private class DoneButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
	        int count = 0;
	        for(int i=0; i<checkBox.length; i++) {
	            if(checkBox[i].isSelected()) count++;
	        }
	        Object[] objects = new Object[count];
	        count = 0;
	        for(int i=0; i<checkBox.length; i++) {
	            if(checkBox[i].isSelected()) objects[count++] = checkBox[i].object;
	        }
	        editor.setValue(objects);
	        editor.firePropertyChange();
	    }
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}