package etomica.virial.GUI.containers;

import java.awt.CardLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.MatteBorder;
import javax.swing.border.SoftBevelBorder;

import etomica.virial.GUI.models.SuperModel;


public class Species extends JPanel implements SubPanelsInterface{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JRadioButton Pure, Binary;
	private JPanel Cards;
	private SingleSpeciesParam OneSpecies;
	private BinaryMixtureParam TwoSpecies;
	
	

	public Species(){
		super();
		initComponents();
	}
	
	
	

	public JRadioButton getPure() {
		return Pure;
	}




	public void setPure(JRadioButton pure) {
		Pure = pure;
	}




	public JRadioButton getBinary() {
		return Binary;
	}




	public void setBinary(JRadioButton binary) {
		Binary = binary;
	}




	public void initComponents(){
		//this.setBackground(Color.yellow);
		this.setLayout(new FlowLayout());
		JPanel ButtonGroup = new JPanel();
		
		this.Cards = new JPanel(new CardLayout());
		ButtonGroup MixtureType = new ButtonGroup();
		//ButtonGroup.setBorder(new SoftBevelBorder(SoftBevelBorder.RAISED));
		this.Pure = new JRadioButton("Pure Species");
		this.Pure.setToolTipText("Select this for pure component");
		this.Binary = new JRadioButton("Binary Mixture");
		MixtureType.add(Pure);
		MixtureType.add(Binary);
		ButtonGroup.add(Pure);
		ButtonGroup.add(Binary);
		//this.add(Pure);
		//this.add(Binary);
		this.add(ButtonGroup);
		this.add(Cards);
		//OneSpecies = new SingleSpeciesParam();
		//TwoSpecies = new BinaryMixtureParam();
		//Cards.add("PURE", OneSpecies);
		//Cards.add("BINARY",TwoSpecies);
		
	}

	public void addChoosePureRadio(ActionListener P){
		this.Pure.addActionListener(P);
	}
	
	public JPanel getCards() {
		return Cards;
	}
	
	public void addCards(SingleSpeciesParam onespecies, BinaryMixtureParam twospecies){
		this.Cards.add("PURE", onespecies);
		this.Cards.add("BINARY", twospecies);
	}




	public void setCards(JPanel cards) {
		Cards = cards;
	}




	public SingleSpeciesParam getOneSpecies() {
		return OneSpecies;
	}




	public void setOneSpecies(SingleSpeciesParam oneSpecies) {
		OneSpecies = oneSpecies;
	}




	public BinaryMixtureParam getTwoSpecies() {
		return TwoSpecies;
	}




	public void setTwoSpecies(BinaryMixtureParam twoSpecies) {
		TwoSpecies = twoSpecies;
	}




	public void addChooseBinaryRadio(ActionListener B){
		this.Binary.addActionListener(B);
	}
	
	
	
}
