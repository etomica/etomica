package etomica;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;

/**
 * Panel for use with a Meter[] property editor.  Permits selection of
 * set of meters from those that are presently instantiated.  Presents
 * an array of checkboxes for the selection.
 */
public class McMoveEditorPanel extends JPanel {

    PropertyEditorSupport editor;
    MyCheckBox[] checkBox; //array of checkboxes with associated mcmove objects
    private static int IDnumber = 0;
    
    public McMoveEditorPanel(PropertyEditorSupport ed) {
        editor = ed;
        
        //discover instantiable MCMove classes
        int counter = 0;
        Class[] objects = etomica.gui.IntrospectionArrays.introspect("MCMove",false);
  //      Object[] objects = Simulation.instance.meterList().toArray();
        
        //construct panel of checkboxes
        JPanel objectPanel = new JPanel(new GridLayout(0,2));
        checkBox = new MyCheckBox[objects.length];
        for(int i=0; i<objects.length; i++) {
            checkBox[i] = new MyCheckBox(objects[i]);
            objectPanel.add(checkBox[i]);
        }
        JButton button = new JButton("Add");
        button.addActionListener(new DoneButtonListener());
                
        add(button);
        add(objectPanel);
    }
    
    /**
     * Simple extension of JCheckBox to permit it to be associated with an object
     */
    private class MyCheckBox extends JCheckBox {
        Class object;
        MyCheckBox(Class obj) {
            super(obj.getName().substring(14));
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
	        MCMove[] moves = new MCMove[count];
	        count = 0;
	        MCMove newMove;
	        for(int i=0; i<checkBox.length; i++) {
	            if(checkBox[i].isSelected()) {
	                try {
	                    newMove = (MCMove)checkBox[i].object.newInstance();
	                    newMove.setName(newMove.getClass().getName().substring(8) + Integer.toString(IDnumber++));
	                } 
	                catch(InstantiationException e) {continue;}
	                catch(IllegalAccessException e) {continue;}
	                moves[count++] = newMove;
	                checkBox[i].setSelected(false);
	            }
	        }
	        editor.setValue(moves); //editor fires property change as part of this call
	    }
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}