package simulate;
import simulate.*;
import simulate.gui.PropertySheet;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.*;

/**
 * Panel for use with a Configuration property editor.
 */
public class ConfigurationEditorPanel extends JPanel {

    PropertyEditorSupport editor;
    MyCheckBox[] checkBox; //array of checkboxes with associated configuration objects
    
    public ConfigurationEditorPanel(PropertyEditorSupport ed) {
        editor = ed;
        int splitPaneHeight = 580;
        int splitPaneWidth = 580;
        setSize(splitPaneWidth, splitPaneHeight);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.anchor = gbc.WEST;
        GridBagLayout gbl = new GridBagLayout();
        javax.swing.JSplitPane splitPane = new javax.swing.JSplitPane();
        add(splitPane);
    /**
     * Lists all of the simulation components corresponding to the respective tabs name.  These are listed
     * as radio buttons.
     */
     javax.swing.JPanel leftPanePanel = new javax.swing.JPanel();
    /**
     * Scrollable pane that holds the leftPanePanel from above
     */
     javax.swing.JScrollPane leftPane = new javax.swing.JScrollPane(leftPanePanel);
     splitPane.setLeftComponent(leftPane);

     PropertySheet rightPanePanel = new PropertySheet(null,0,0);
     javax.swing.JScrollPane rightPane = new javax.swing.JScrollPane(rightPanePanel);
     splitPane.setRightComponent(rightPane);
     
        int leftPanelHeight = splitPaneHeight;
        int leftPanelWidth = (int)(0.4*splitPaneWidth)-10;//10 accounts for width of slider bar
        leftPanePanel.setLayout(gbl);
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
    
        //discover configuration presently instantiated
/*        int counter = 0;
        Object[] objects = new Object[
            SimEditorTabMenu.getConfigurationEditor().getComponentList().size()];
        for (int i = 0; i < SimEditorTabMenu.getConfigurationEditor().getComponentList().size(); i++)
            objects[counter++] = SimEditorTabMenu.getConfigurationEditor().getComponentList().elementAt(i);
        
        //construct panel of checkboxes
        Configuration currentConfiguration = (Configuration)editor.getValue();
        JPanel objectPanel = new JPanel(new GridLayout(0,2));
        checkBox = new MyCheckBox[objects.length];
        for(int i=0; i<objects.length; i++) {
            checkBox[i] = new MyCheckBox(objects[i]);
            objectPanel.add(checkBox[i]);
            for(int j=0; j<currentConfiguration.length; j++) {  //put check in box if already in list of configuration
                if(currentConfiguration[j] == objects[i]) checkBox[i].setSelected(true);
            }
        }
        JButton button = new JButton("OK");
        button.addActionListener(new DoneButtonListener());
                
        add(button);
        add(objectPanel);
        */
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
//	        Configuration configuration = new Configuration[count];
	        count = 0;
	        for(int i=0; i<checkBox.length; i++) {
//	            if(checkBox[i].isSelected()) configuration[count++] = (Configuration)checkBox[i].object;
	        }
//	        editor.setValue(configuration);
	        editor.firePropertyChange();
	    }
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
}