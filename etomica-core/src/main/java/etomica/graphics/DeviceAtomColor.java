/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.graphics;

import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;

public class DeviceAtomColor extends Device implements ChangeListener {

    protected final JPanel panel;
    protected final JColorChooser jcc;
    protected final AtomType type;
    protected final ColorSchemeByType colorScheme;
    protected final IAction repaintAction;

    public DeviceAtomColor(Controller controller, ColorSchemeByType colorScheme, AtomType type, IAction repaintAction) {
        super(controller);
        this.type = type;
        this.colorScheme = colorScheme;
        panel = new JPanel(new java.awt.BorderLayout());

        // Set up color chooser for setting text color
        jcc = new JColorChooser(Color.RED);
        jcc.getSelectionModel().addChangeListener(this);
        jcc.setBorder(BorderFactory.createTitledBorder(
                "Choose Atom Color"));
        jcc.setPreviewPanel(new JPanel());
        panel.add(jcc, java.awt.BorderLayout.CENTER);
        this.repaintAction = repaintAction;
    }

    @Override
    public Component graphic() {
        return panel;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        Color newColor = jcc.getColor();
        colorScheme.setColor(type, newColor);
        if (repaintAction!=null) repaintAction.actionPerformed();
    }
}
