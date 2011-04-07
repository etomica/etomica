package etomica.virial.GUI.containers;

import java.awt.Color;

import javax.swing.JPanel;

import etomica.virial.GUI.models.SuperModel;


public class MoveParam extends JPanel implements SubPanelsInterface{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	MoveParam(){
		super();
		this.initComponents();
	}
	
	
	public void initComponents(){
		this.setBackground(Color.blue);
	}

}
