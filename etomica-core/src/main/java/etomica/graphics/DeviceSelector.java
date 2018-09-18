/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.HashMap;

import etomica.action.IAction;
import etomica.action.activity.Controller;

/**
 * Generic Device that has a combo box.  Selecting an item from the combo box
 * will invoke an action.  Each item has its own action.
 */
public class DeviceSelector extends Device {

     public DeviceSelector(Controller controller) {
        super(controller);
        selector = new javax.swing.JComboBox();
        label = new javax.swing.JLabel("");

        panel = new javax.swing.JPanel(new java.awt.BorderLayout(0,1));
        panel.add(label, java.awt.BorderLayout.NORTH);
        panel.add(selector, java.awt.BorderLayout.SOUTH);
        panel.setBorder(new javax.swing.border.EmptyBorder(3,3,3,3));
        actionHash = new HashMap<String,IAction>();
        
        //listener to combo box gets value and initiates action
        selector.addItemListener( new ItemListener() {
            public void itemStateChanged(ItemEvent event) {
                if (event.getStateChange() == ItemEvent.DESELECTED) return;
                IAction action = actionHash.get(event.getItem());
                doAction(action);
            }
        });
    }

    public javax.swing.JComboBox getSelector() {return selector;}
    public javax.swing.JLabel getLabel() {return label;}
    
    public void setLabel(String newLabel) {
        label.setText(newLabel);
    }

    public void addOption(String actionText, IAction action) {
        selector.addItem(actionText);
        actionHash.put(actionText, action);
    }

    public void removeOption(String actionText) {
        selector.removeItem(actionText);
        actionHash.remove(actionText);
    }

    /**
     * Sets the i-th item (counting from 0) in the list as the one currently selected.
     * Argument outside bounds will result in an exception.
     */
    public void setSelected(int i) {
        if(i < 0 || i >= selector.getItemCount()) throw new IllegalArgumentException("Out of range");
        selector.setSelectedIndex(i);
    }
    
    /**
     * Returns the GUI element for display in the simulation.
     * Consists of a combo box used for the selector.
     */
    public java.awt.Component graphic(Object obj) {
        return panel;
    }
    
    private javax.swing.JComboBox selector;
    private javax.swing.JLabel label;
    private javax.swing.JPanel panel;
    protected final HashMap<String,IAction> actionHash;
}
