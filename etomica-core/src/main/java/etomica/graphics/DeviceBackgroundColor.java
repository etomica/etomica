/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.graphics;

import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Device to set background color of configuration display.
 *
 * Has been tested only for 2D simulation; uncertain if it will work for 3D
 */
public class DeviceBackgroundColor extends Device implements ChangeListener {

    protected final JPanel panel;
    protected final JColorChooser jcc;
    protected final DisplayCanvas displayCanvas;
    protected final IAction repaintAction;

    public DeviceBackgroundColor(Controller controller, SimulationGraphic sim, Box box, IAction repaintAction) {
        super(controller);
        this.displayCanvas = sim.getDisplayBox(box).canvas;
        panel = new JPanel(new BorderLayout());

        jcc = new JColorChooser(Color.WHITE);
        jcc.getSelectionModel().addChangeListener(this);
        jcc.setBorder(BorderFactory.createTitledBorder(
                "Choose Background Color"));
        jcc.setPreviewPanel(new JPanel());
        panel.add(jcc, BorderLayout.CENTER);
        this.repaintAction = repaintAction;
    }

    @Override
    public Component graphic() {
        return panel;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        Color newColor = jcc.getColor();
        displayCanvas.setBackground(newColor);
        if (repaintAction!=null) repaintAction.actionPerformed();
    }

    public class Button extends JButton {

        public Button() {
            super("Background");

            addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    JFrame f = new JFrame();
                    f.getContentPane().add(DeviceBackgroundColor.this.graphic());
                    f.pack();
                    f.setTitle("Background Color");
                    f.setVisible(true);
                }
            });
        }
    }
}
