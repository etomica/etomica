
// Support for PropertyEditors that use tags.

package etomica.gui;

import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

class PropertySelector extends JComboBox implements ItemListener {

    PropertySelector(PropertyEditor pe) {
	    editor = pe;
	    String tags[] = editor.getTags();
	    String current = editor.getAsText();
	    int iCurrent = 0;
	    DefaultComboBoxModel componentList = new DefaultComboBoxModel();
	    ((JComboBox)this).setModel(componentList);
	    for (int i = 0; i < tags.length; i++) {
	        ((JComboBox)this).addItem(tags[i]);
	        if(current == tags[i]) iCurrent = i;
	    }
	    ((JComboBox)this).setSelectedIndex(iCurrent);
    	addItemListener(this);
    	setBorder(PropertySheet.EMPTY_BORDER);
    }

    public void itemStateChanged(ItemEvent evt) {
	    editor.setAsText(((String)((JComboBox)evt.getSource()).getSelectedItem()));
	    transferFocus();
    }

    transient PropertyEditor editor;    
}
