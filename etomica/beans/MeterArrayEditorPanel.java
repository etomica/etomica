package simulate;
import simulate.*;
import simulate.gui.SimEditorTabMenu;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;

/**
 * Panel for use with a Meter[] property editor.  Permits selection of
 * set of meters from those that are presently instantiated.  Presents
 * an array of checkboxes for the selection.
 */
public class MeterArrayEditorPanel extends JPanel {

    PropertyEditorSupport editor;
    MyCheckBox[] checkBox; //array of checkboxes with associated meter objects
    
    public MeterArrayEditorPanel(PropertyEditorSupport ed) {
        editor = ed;
        
        //discover meters presently instantiated
        int counter = 0;
        Object[] objects = new Object[
            SimEditorTabMenu.getMeterEditor().getComponentList().size()];
        for (int i = 0; i < SimEditorTabMenu.getMeterEditor().getComponentList().size(); i++)
            objects[counter++] = SimEditorTabMenu.getMeterEditor().getComponentList().elementAt(i);
        
        //construct panel of checkboxes
        Meter[] currentMeters = (Meter[])editor.getValue();
        JPanel objectPanel = new JPanel(new GridLayout(0,2));
        checkBox = new MyCheckBox[objects.length];
        for(int i=0; i<objects.length; i++) {
            checkBox[i] = new MyCheckBox(objects[i]);
            objectPanel.add(checkBox[i]);
            for(int j=0; j<currentMeters.length; j++) {  //put check in box if already in list of meters
                if(currentMeters[j] == objects[i]) checkBox[i].setSelected(true);
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
	        Meter[] meters = new Meter[count];
	        count = 0;
	        for(int i=0; i<checkBox.length; i++) {
	            if(checkBox[i].isSelected()) meters[count++] = (Meter)checkBox[i].object;
	        }
	        editor.setValue(meters);
	        editor.firePropertyChange();
	    }
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}