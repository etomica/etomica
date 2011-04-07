package etomica.virial.GUI.containers;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.SoftBevelBorder;

public class Species1 implements ActionListener{
	
	private JButton AddSpecies;
	private JFrame frame;
	private String[] SpeciesList = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane","Naphthalene"};
	private String[] IntialList = {"---No species selected---","Press \'Add\' to select a species"};
	private JList Specieslist;
	private JScrollPane listScroller;
	private JTextArea InitialSetup;
	
	
	
	
	Species1(){
	frame = new JFrame();
	
	JPanel panel = new JPanel();
	panel.setLayout(new FlowLayout());
	frame.add(panel);
	
	AddSpecies = new JButton("Add Species");
	AddSpecies.addActionListener(this);
	
	panel.add(AddSpecies);
	
	Specieslist = new JList(IntialList);
	Specieslist.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
	Specieslist.setLayoutOrientation(JList.HORIZONTAL_WRAP);
	Specieslist.setVisibleRowCount(-1);
	
 
    listScroller = new JScrollPane(Specieslist);
    listScroller.setPreferredSize(new Dimension(300,100));
    listScroller.setAlignmentX(JComponent.LEFT_ALIGNMENT);
     
    JPanel listPane = new JPanel();
    listPane.setLayout(new BoxLayout(listPane, BoxLayout.PAGE_AXIS));
    listPane.setBorder(new SoftBevelBorder(SoftBevelBorder.RAISED));
    JLabel label = new JLabel("Species");
    label.setLabelFor(Specieslist);
    listPane.add(label);
    listPane.add(Box.createRigidArea(new Dimension(0,5)));
    listPane.add(listScroller);
    listPane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
     
    panel.add(listPane);
	
	}

	public void actionPerformed(ActionEvent e){
		
		Specieslist.setListData(SpeciesList);
		
	}

	public JList getList() {
		return Specieslist;
	}

	public void setList(JList list) {
		this.Specieslist = list;
	}

	public JButton getAddSpecies() {
		return AddSpecies;
	}


	public void setAddSpecies(JButton addSpecies) {
		AddSpecies = addSpecies;
	}
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				Species1 s = new Species1();
				s.frame.setVisible(true);
				s.frame.setResizable(true);
				
				
				
			}
		});
    }

}
