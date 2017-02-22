/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.containers;

import javax.swing.JPanel;

import java.awt.Dimension;

public class PanelMainFrame extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public PanelMainFrame(int Width, int Height){
		super();
		this.setAlignmentX((float) 0.5);
		this.setAlignmentY((float) 0.5);
		this.setPreferredSize(new Dimension(Width, Height));
	}
}
