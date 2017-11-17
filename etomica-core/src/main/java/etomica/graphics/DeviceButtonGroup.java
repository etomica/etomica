/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.action.activity.Controller;

/**
 * Add etomica device for a group of radio buttons.  Each button is associated
 * with an IAction.  The buttons can be referenced by their label; their
 * selection status can be set or queried.
 *
 * @author Andrew Schultz
 */
public class DeviceButtonGroup extends Device {

    protected final ButtonGroup buttonGroup;
    protected final ArrayList<JRadioButton> buttons;
    protected final ArrayList<IAction> buttonActions;
    protected final HashMap<String,Integer> labelHash;
    protected final JPanel buttonPanel;
    protected final int cols;
    protected IAction targetAction;
    
    public DeviceButtonGroup(Controller controller, int cols) {
        super(controller);
        buttonGroup = new ButtonGroup();
        buttons = new ArrayList<JRadioButton>();
        buttonActions = new ArrayList<IAction>();
        labelHash = new HashMap<String, Integer>();
        this.cols = cols;
        
        buttonPanel = new JPanel(new GridBagLayout());
    }

    /**
     * Sets the label for the group of radio buttons
     */
    public void setLabel(String label) {
        buttonPanel.setBorder(new TitledBorder(null, label, TitledBorder.CENTER, TitledBorder.TOP));
    }

    /**
     * Adds a radio button with the given label, which (when selected) performs
     * the given action.
     */
    public JRadioButton addButton(String label, final IAction action) {
        int id = buttons.size();
        JRadioButton newButton = new JRadioButton(label);
        buttons.add(newButton);
        buttonGroup.add(newButton);
        labelHash.put(label, id);
        
        GridBagConstraints gbc1 = new GridBagConstraints();
        gbc1.gridx = id % cols;  gbc1.gridy = id/cols;
        gbc1.gridwidth = 1;
        buttonPanel.add(newButton, gbc1);
        if (action != null) {
            ActionListener buttonAction = new ActionListener() {
                public void actionPerformed(ActionEvent evt) {
                    try {
                        doAction(action);
                    }
                    catch (RuntimeException e) {
                        System.err.println(e+" "+e.getMessage());
                    }
                }
            };
            newButton.addActionListener(buttonAction);
            while (id > buttonActions.size()+1) {
                buttonActions.add(null);
            }
            buttonActions.add(action);
        }
        return newButton;
    }

    /**
     * Returns true if the button of the given label is selected.
     */
    public boolean isSelected(String label) {
        int i = labelHash.get(label);
        return buttons.get(i).isSelected();
    }
    
    /**
     * Selects the button with the given label.
     */
    public void setSelected(String label) {
        setSelected(labelHash.get(label));
    }
    
    /**
     * Selects ith radio button.
     */
    protected void setSelected(int i) {
        buttons.get(i).setSelected(true);
        IAction action = buttonActions.get(i);
        if (action == null) return;
        try {
            doAction(action);
        }
        catch (RuntimeException e) {
            System.err.println(e+" "+e.getMessage());
        }
    }
    
    public Component graphic(Object obj) {
        return buttonPanel;
    }

}
