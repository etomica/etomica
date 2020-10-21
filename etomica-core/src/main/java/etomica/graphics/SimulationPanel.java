/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

/**
 * Main panel used to display the graphical elements of a simulation.  
 * An instance of this panel is first constructed when the panel() method
 * of a Simulation instance is invoked.  This method is called when the panel
 * is added to an applet or a window, and it is also invoked by the default
 * Mediator class when a Display or Device are added to a Simulation.
 *
 * @author David Kofke
 */

/* History
 * 
 * 01/03/03 (DAK) added constructor that omits addition of internal panels
 */
 
public class SimulationPanel extends JPanel {
    
	private static final GridBagConstraints vertGBC = new GridBagConstraints();
	private static final GridBagConstraints horizGBC = new GridBagConstraints();

	public final JPanel toolbarPanel = new JPanel();
	public final JPanel controlPanel = new JPanel();
    public final JPanel graphicsPanel = new JPanel();
    public final JPanel footerPanel = new JPanel();
    public final JPanel plotPanel = new JPanel();
    public final JPanel metricPanel = new JPanel();
    public final JTabbedPane tabbedPane = new JTabbedPane();
    public final DefaultToolbar toolbar;

    public SimulationPanel() {
    	this("");
    }

    public SimulationPanel(String appName) {
        toolbarPanel.setLayout(new GridLayout());
        controlPanel.setLayout(new GridBagLayout());
        footerPanel.setLayout(new GridBagLayout());
        graphicsPanel.setLayout(new BorderLayout(5, 5));
        plotPanel.setLayout(new GridBagLayout());
    	metricPanel.setLayout(new GridBagLayout());
        setLayout(new BorderLayout(10, 10));
        toolbar = new DefaultToolbar(this, appName);
        toolbarPanel.add(toolbar.graphic());
        add(toolbarPanel, BorderLayout.NORTH);
        add(controlPanel, BorderLayout.WEST);
        add(graphicsPanel);
        add(footerPanel, BorderLayout.SOUTH);
        add(plotPanel, BorderLayout.EAST);
    }


    public static GridBagConstraints getVertGBC() {
    	vertGBC.gridx = 0;
    	vertGBC.gridy = GridBagConstraints.RELATIVE;
        vertGBC.insets = new java.awt.Insets(3,1,3,1);
    	return vertGBC;
    }
    
    public static GridBagConstraints getHorizGBC() {
    	horizGBC.gridx = GridBagConstraints.RELATIVE;
    	horizGBC.gridy = 0;
        horizGBC.insets = new java.awt.Insets(3,1,3,1);
    	return horizGBC;
    }

}//end of SimulationPanel
