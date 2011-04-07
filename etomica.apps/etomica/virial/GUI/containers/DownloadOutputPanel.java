package etomica.virial.GUI.containers;

import java.awt.BorderLayout;
import java.awt.Color;


import javax.swing.BorderFactory;
import javax.swing.JPanel;

public class DownloadOutputPanel extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public DownloadOutputPanel(String PanelTitle){
		super(new BorderLayout());
		this.setBackground(Color.yellow);
		this.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder(PanelTitle),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	}

	
	
}
