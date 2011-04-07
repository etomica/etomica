package etomica.virial.GUI.containers;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JFrame;



public class MainFrame extends JFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	public MainFrame(String FrameLabel){
		super(FrameLabel);
		
	}
	
	

	public void SetFrameProperties(){
		this.setVisible(true);
		this.setResizable(true);
		
		//this.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);  
		this.pack();
	}


	
}
