package etomica.virial.GUI.containers;

import java.awt.FlowLayout;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class DialogBoxPanel extends JPanel{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JLabel Message;
	private JButton CloseWindow;

	DialogBoxPanel(){
		Message = new JLabel("This Simulation cannot be run now. Please try again later");
		CloseWindow = new JButton("Ok");
		this.setLayout(new FlowLayout());
		this.add(Message);
		this.add(CloseWindow);

	}
	
	public void addChooseOkayButtonListener(ActionListener Comp) {
		CloseWindow.addActionListener(Comp);
	}

}
