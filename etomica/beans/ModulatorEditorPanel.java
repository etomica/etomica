package simulate;
import simulate.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;
//import com.symantec.itools.vcafe.openapi.*;
import simulate.utility.StringUtility;

public class ModulatorEditorPanel extends JPanel {

    PropertyEditorSupport editor;
    JPanel propertyPanel;
///    String currentObject;
    Object currentObject;
    
    public ModulatorEditorPanel(PropertyEditorSupport ed) {
        editor = ed;
        int counter = 0;
        
        Object[] objects = Simulation.instance.allElements().toArray();

        setLayout(new GridLayout(1,2));//this panel is made of two subpanels side-by-side
        
        ButtonGroup buttonGroup = new ButtonGroup();
        JPanel objectPanel = new JPanel(new GridLayout(0,1));
        for(int i=0; i<objects.length; i++) {
            MyRadioButton objectButton = new MyRadioButton(objects[i],i==0);
            buttonGroup.add(objectButton);
            objectPanel.add(objectButton);
            objectButton.addActionListener(new ObjectButtonListener());
        }
        add(objectPanel);
        
        propertyPanel = new JPanel(new GridLayout(0,1));//number of rows to be determined, one column
        add(propertyPanel);
    }
    
    /**
     * Simple extension of JRadioButton to permit it to be associated with an object
     */
    private class MyRadioButton extends JRadioButton {
        Object object;
        MyRadioButton(Object obj, boolean selected) {
            super(obj.toString(),selected);
            object = obj;
        }
    }
    
    /**
     * ActionListener class called when an object is selected.
     * Action is to set up property list for object.
     */
	private class ObjectButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
            BeanInfo bi = null;
	        
	        
	        Object object = (Object)((MyRadioButton)evt.getSource()).object;
 ///           currentObject = object.toString();
            currentObject = object;
	        
	        try {
	            bi = Introspector.getBeanInfo(object.getClass());
	        } 
	        catch (IntrospectionException ex) {
	            System.out.println("Couldn't introspect");
	        }
	        
	        //VisualProperty[] properties = object.getProperties();
            PropertyDescriptor[] properties = bi.getPropertyDescriptors();
            ButtonGroup buttonGroup = new ButtonGroup();
            propertyPanel.removeAll();
            for(int i=0; i<properties.length; i++) {
                if(properties[i].getName().length() <= 1) continue;
                JRadioButton propertyButton = new JRadioButton(StringUtility.decapitalize(properties[i].getName()),i==0);
                buttonGroup.add(propertyButton);
                propertyPanel.add(propertyButton);
                propertyButton.addActionListener(new PropertyButtonListener());
            }
            validate();
            propertyPanel.repaint();
	    }
	}
	
    /**
     * ActionListener class called when a property is selected.
     * Action is to set value in editor and fire property change event.
     */
	private class PropertyButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
/*	        editor.setValue(new String[] 
	            {
	                currentObject,
	                ((JRadioButton)evt.getSource()).getText()
	            });
	            */
	        editor.setValue(new Modulator(currentObject, ((JRadioButton)evt.getSource()).getText()));
	        editor.firePropertyChange();
	    }
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}