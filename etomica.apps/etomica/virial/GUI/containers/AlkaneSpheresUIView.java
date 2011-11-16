package etomica.virial.GUI.containers;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import etomica.space1d.Space1D;

public class AlkaneSpheresUIView  extends JPanel {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public JButton CloseWindow;
	public JButton SaveValues;
	
	public JTextField NoOfSpheres;
	//private static final NAlkaneSpheresParameter INSTANCE = new NAlkaneSpheresParameter();
	
	public AlkaneSpheresUIView(){
		
		
		CloseWindow = new JButton("Cancel");
		SaveValues = new JButton("Save");
		GridBagLayout gridbagOtherParam = new GridBagLayout();
		this.setLayout(gridbagOtherParam);
	    
		NoOfSpheres = new JTextField();
		NoOfSpheres.setColumns(10);
		NoOfSpheres.setSize(new Dimension(5,20));
		NoOfSpheres.setText("4");
		JLabel NLabel = new JLabel("Cn (where n is 4 or greater)");
		NLabel.setLabelFor(NoOfSpheres);
		
		
		Border compoundField1;
		compoundField1 = BorderFactory.createCompoundBorder(
	    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		NoOfSpheres.setBorder(compoundField1);
		
		JComponent[] OtherParamLeft= {NLabel,SaveValues};
		JComponent[] OtherParamRight = { NoOfSpheres,CloseWindow};
	    addLeftRightComponents(OtherParamLeft,OtherParamRight,gridbagOtherParam,this);
		
		
	}
	/*
	public static NAlkaneSpheresParameter getInstance() {
		return INSTANCE;
		
	}*/
	private void addLeftRightComponents(JComponent[] ComponentLeft,
			JComponent[] ComponentRight,
            GridBagLayout gridbag,
            Container container) {
		
			GridBagConstraints c = new GridBagConstraints();
			c.anchor = GridBagConstraints.NORTHEAST;
			int numLabels = ComponentLeft.length;

			for (int i = 0; i < numLabels; i++) {
				c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
				c.fill = GridBagConstraints.NONE;      //reset to default
				c.weightx = 0.0;   
				c.insets = new Insets(0,10,0,0);//reset to default
				container.add(ComponentLeft[i], c);

				c.gridwidth = GridBagConstraints.REMAINDER;     //end row
				c.fill = GridBagConstraints.HORIZONTAL;
				c.weightx = 0.0;
				c.insets = new Insets(0,10,3,0);
				container.add(ComponentRight[i], c);
			}
	}
	public JTextField getNoOfSpheres() {
		return NoOfSpheres;
	}

	public void setNoOfSpheres(JTextField noOfSpheres) {
		NoOfSpheres = noOfSpheres;
	}

	public JButton getCloseWindow() {
		return CloseWindow;
	}

	public void setCloseWindow(JButton closeWindow) {
		CloseWindow = closeWindow;
	}

	public JButton getSaveValues() {
		return SaveValues;
	}

	public void setSaveValues(JButton saveValues) {
		SaveValues = saveValues;
	}
	
}
