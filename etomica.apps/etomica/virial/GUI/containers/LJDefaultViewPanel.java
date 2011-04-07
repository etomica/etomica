package etomica.virial.GUI.containers;

import java.awt.CardLayout;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.text.DecimalFormat;

import javax.swing.JButton;
import javax.swing.JComponent;

import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.PlainDocument;

public class LJDefaultViewPanel extends JPanel{
	
	private JLabel Sigma1;
	private JLabel Sigma2;
	private JLabel Epsilon1;
	private JLabel Epsilon2;
	private JLabel Moment;
	private JLabel BondLength;
	
	private JLabel SigmaAA;
	private JLabel SigmaBB;
	private JLabel SigmaAB;
	
	private JLabel EpsilonAA;
	private JLabel EpsilonBB;
	private JLabel EpsilonAB;
	
	private JTextField SigmaField1;
	private JTextField SigmaField2;
	private JTextField SigmaFieldAA;
	private JTextField SigmaFieldBB;
	private JTextField SigmaFieldAB;
	
	private JTextField EpsilonField1;
	private JTextField EpsilonField2;
	private JTextField EpsilonFieldAA;
	private JTextField EpsilonFieldAB;
	private JTextField EpsilonFieldBB;
	
	private JTextField MomentField;
	private JTextField BondLField;
	
	private JPanel LJDefaultViewCards;
	private JButton Help1;
	private JButton Help2;
	private JButton Help3;
	private JButton Help4;
	
	
	
	public LJDefaultViewPanel(){
		super();
		
		LJDefaultViewCards = new JPanel(new CardLayout());
		
		Sigma1 = new JLabel("Sigma");
		Epsilon1 = new JLabel("Epsilon");
		
		SigmaField1 = new JTextField(10);
		EpsilonField1 = new JTextField(10);
		
		
		
		Sigma2 = new JLabel("Sigma");
		SigmaField2 = new JTextField(10);
		
		Epsilon2 = new JLabel("Epsilon"); 
		EpsilonField2 = new JTextField(10);
		Moment = new JLabel("Moment");
		BondLength = new JLabel("Bond Length"); 
		MomentField = new JTextField(10);
		BondLField = new JTextField(10);
		
		
		SigmaAA = new JLabel("Sigma11");
		SigmaBB = new JLabel("Sigma22");
		SigmaAB = new JLabel("Sigma12");
		
		EpsilonAA = new JLabel("Epsilon11");
		EpsilonBB = new JLabel("Epsilon22");
		EpsilonAB = new JLabel("Epsilon12");
		
		
		SigmaFieldAA = new JTextField(10);
		SigmaFieldAB = new JTextField(10);
		SigmaFieldBB = new JTextField(10);
		
	
		
		
		EpsilonFieldAA = new JTextField(10);
		EpsilonFieldBB = new JTextField(10);
		EpsilonFieldAB = new JTextField(10);
		
		
		
		
		
		JPanel LJ1 = new JPanel();
		GridBagLayout gridbag1 = new GridBagLayout();
		LJ1.setLayout(gridbag1);
		
		JLabel[] labelsLJ1 = {Sigma1,Epsilon1};
		JComponent[] containerLJ1 = {SigmaField1,EpsilonField1};
		addLabelWithComponents(labelsLJ1,containerLJ1,gridbag1,LJ1);
		
		
		
		
		JPanel LJ2 = new JPanel();
		GridBagLayout gridbag2 = new GridBagLayout();
		LJ2.setLayout(gridbag2);
		
		JLabel[] labelsLJ2 = {Sigma2,Epsilon2,Moment,BondLength};
		JComponent[] containerLJ2 = {SigmaField2,EpsilonField2,MomentField,BondLField};
		addLabelWithComponents(labelsLJ2,containerLJ2,gridbag2,LJ2);
		
		
		
		JPanel LJMix = new JPanel();
		GridBagLayout gridbagMix = new GridBagLayout();
		LJMix.setLayout(gridbagMix);
		
		JLabel[] labelsLJMix = {SigmaAA,SigmaAB, SigmaBB,EpsilonAA,EpsilonAB,EpsilonBB};
		JComponent[] containerLJMix  = {SigmaFieldAA,SigmaFieldAB,SigmaFieldBB,EpsilonFieldAA,EpsilonFieldAB,EpsilonFieldBB};
		addLabelWithComponents(labelsLJMix,containerLJMix,gridbagMix,LJMix);
		
		LJDefaultViewCards.add("LJ1",LJ1);
		LJDefaultViewCards.add("LJ2",LJ2);
		LJDefaultViewCards.add("LJMix",LJMix);
		this.add(LJDefaultViewCards);
		
		
	}
	
	
	
	
	public JTextField getSigmaField1() {
		return SigmaField1;
	}

	public void setSigmaField1(JTextField sigmaField) {
		SigmaField1 = sigmaField;
	}
	
	public JTextField getSigmaField2() {
		return SigmaField2;
	}

	public void setSigmaField2(JTextField sigmaField) {
		SigmaField2 = sigmaField;
	}

	public JTextField getSigmaFieldAA() {
		return SigmaFieldAA;
	}

	public void setSigmaFieldAA(JTextField sigmaFieldAA) {
		SigmaFieldAA = sigmaFieldAA;
	}

	public JTextField getSigmaFieldBB() {
		return SigmaFieldBB;
	}

	public void setSigmaFieldBB(JTextField sigmaFieldBB) {
		SigmaFieldBB = sigmaFieldBB;
	}

	public JTextField getSigmaFieldAB() {
		return SigmaFieldAB;
	}

	public void setSigmaFieldAB(JTextField sigmaFieldAB) {
		SigmaFieldAB = sigmaFieldAB;
	}

	public JTextField getEpsilonField1() {
		return EpsilonField1;
	}

	public void setEpsilonField1(JTextField epsilonField1) {
		EpsilonField1 = epsilonField1;
	}
	
	public JTextField getEpsilonField2() {
		return EpsilonField2;
	}

	public void setEpsilonField2(JTextField epsilonField) {
		EpsilonField2 = epsilonField;
	}

	public JTextField getEpsilonFieldAA() {
		return EpsilonFieldAA;
	}

	public void setEpsilonFieldAA(JTextField epsilonFieldAA) {
		EpsilonFieldAA = epsilonFieldAA;
	}

	public JTextField getEpsilonFieldAB() {
		return EpsilonFieldAB;
	}

	public void setEpsilonFieldAB(JTextField epsilonFieldAB) {
		EpsilonFieldAB = epsilonFieldAB;
	}

	public JTextField getEpsilonFieldBB() {
		return EpsilonFieldBB;
	}

	public void setEpsilonFieldBB(JTextField epsilonFieldBB) {
		EpsilonFieldBB = epsilonFieldBB;
	}

	public JTextField getMomentField() {
		return MomentField;
	}

	public void setMomentField(JTextField momentField) {
		MomentField = momentField;
	}

	public JTextField getBondLField() {
		return BondLField;
	}

	public void setBondLField(JTextField bondLField) {
		BondLField = bondLField;
	}

	public JPanel getLJDefaultViewCards() {
		return LJDefaultViewCards;
	}

	public void setLJDefaultViewCards(JPanel lJDefaultViewCards) {
		LJDefaultViewCards = lJDefaultViewCards;
	}

	private void addLabelWithComponents(JLabel[] labels,
            JComponent[] components,
            GridBagLayout gridbag,
            Container container) {

			GridBagConstraints c = new GridBagConstraints();
			//c.anchor = GridBagConstraints.WEST;
			
			
			int numLabels = labels.length;

			for (int i = 0; i < numLabels; i++) {
		
			c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
			//c.fill = GridBagConstraints.NONE;      //reset to default
			c.fill = GridBagConstraints.BOTH;
			c.insets = new Insets(10,10,10,10);
			c.ipady = 5;
			c.anchor = GridBagConstraints.NORTHEAST;
			c.weightx = 0.5;     
			
			
			//reset to default
			container.add(labels[i], c);

			c.gridwidth = GridBagConstraints.REMAINDER;     //end row
			c.fill = GridBagConstraints.BOTH;
			c.anchor = GridBagConstraints.EAST;
			//c.weightx = 1.0;
			container.add(components[i], c);}
		}
	
}
