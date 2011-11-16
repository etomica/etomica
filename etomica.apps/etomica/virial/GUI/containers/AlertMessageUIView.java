package etomica.virial.GUI.containers;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class AlertMessageUIView extends JPanel{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JLabel Message;
	
	private JButton CloseWindow;

	AlertMessageUIView(JFrame frame, String message){
		//JFrame MessageFrame = new JFrame("Alert");
		frame.setMinimumSize(new Dimension(400,100));
		Message = new JLabel(message);
		CloseWindow = new JButton("Ok");
		
		frame.setLayout(new FlowLayout());
		frame.add(Message);
		frame.add(CloseWindow);
		frame.setVisible(true);

	}

	

	public JButton getCloseWindow() {
		return CloseWindow;
	}

	public void setCloseWindow(JButton closeWindow) {
		CloseWindow = closeWindow;
	}
	
	

}
