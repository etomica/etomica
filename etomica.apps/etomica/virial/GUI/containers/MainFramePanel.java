package etomica.virial.GUI.containers;

import javax.swing.JPanel;

import java.awt.Dimension;

public class MainFramePanel extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public MainFramePanel(int Width, int Height){
		super();
		this.setAlignmentX((float) 0.5);
		this.setAlignmentY((float) 0.5);
		this.setPreferredSize(new Dimension(Width, Height));
	}
}
