package etomica.virial.GUI.containers;

import java.awt.BorderLayout;
import java.awt.Color;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;



public class ConsolePanel extends JPanel{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public ConsolePanel(String PanelTitle, JFrame MainFrame){
		super(new BorderLayout());
		this.setBackground(Color.red);
		this.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder(PanelTitle),
                BorderFactory.createEmptyBorder(5,5,5,5)));
		//new ConsoleFrame(this,MainFrame);

	}

}
