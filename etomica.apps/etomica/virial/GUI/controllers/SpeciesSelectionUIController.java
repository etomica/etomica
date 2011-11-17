package etomica.virial.GUI.controllers;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.SwingUtilities;

import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import etomica.virial.GUI.containers.AlertMessageUIView;
import etomica.virial.GUI.containers.AlkaneSpheresUIView;
import etomica.virial.GUI.containers.SimulationEnvironmentUIView;

import etomica.virial.GUI.containers.SpeciesSelectionUIView;
import etomica.virial.GUI.models.AlkaneLengthSelectionDM;
import etomica.virial.GUI.models.SpeciesDM;
import etomica.virial.GUI.models.SpeciesSelectionDM;


public class SpeciesSelectionUIController {
		
	private SpeciesSelectionUIView speciesSelectionView;
	private AlkaneSpheresUIView alkaneLengthSelectionView1;
	private AlkaneSpheresUIView alkaneLengthSelectionView2;
	private AlertMessageUIView alertMessageView;
	private SimulationEnvironmentUIView simEnvView;
	
	
	//DataModels required
	private SpeciesSelectionDM speciesSelectionDM;
	private SpeciesDM speciesDM;
	private AlkaneLengthSelectionDM alkaneLengthSelectionDM;
	
	SpeciesSelectionUIController(SpeciesSelectionUIView speciesSelectionV, 
										SpeciesSelectionDM  speciesSelectionDM, SpeciesDM speciesDM ){
		this.speciesSelectionView = speciesSelectionV;
		this.speciesSelectionDM = speciesSelectionDM;
		this.speciesDM = speciesDM;
		
		speciesSelectionView.addSpeciesButtonListener(new SpeciesAddButtonListener());
		speciesSelectionView.addSpeciesJListListener(new SpeciesJListChangeListener());
		speciesSelectionView.addPotentialButtonListener(new PotentialAddButtonListener());
		speciesSelectionView.addPotentialsJListListener(new PotentialsJListChangeListener());
		
	}
	
	class SpeciesAddButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			speciesSelectionView.getSpeciesJList().setListData(speciesSelectionDM.getSpeciesList());
			speciesSelectionView.getPotentialJList().setEnabled(true);
			speciesSelectionView.getAddPotential().setEnabled(true);
			//This AddFlag is true when a Different Species is to be added to the list
			//When AddFlag is false, we re changing the species selection for potential 1, thats why we give a new instance
			//of the the sList ArrayList
			if(!speciesDM.isSpecies1Added()){
				//sList = new SpeciesList();
				speciesSelectionView.setsList(new ArrayList<String>());
			}
			
		}
	 
	 }
	
	class SpeciesJListChangeListener implements ListSelectionListener
	 {

		@Override
		public void valueChanged(ListSelectionEvent e) {
			// TODO Auto-generated method stub
			if(!speciesSelectionView.getSpeciesJList().getValueIsAdjusting()){
				if(speciesSelectionView.getsList().size() > 1 && !speciesDM.isSpecies1Added()){
					speciesSelectionView.getsList().clear();
					speciesSelectionView.getLJP1().removeData();
					speciesSelectionView.getLJP2().removeData();
				}
				speciesDM.setSpeciesSelectionChanged(false);
				speciesSelectionDM.setSpeciesIdentifier(null);
				speciesSelectionView.getPotentialJList().setListData(speciesSelectionDM.getIntialpotentialdisplaylist());
			}
		}
	 
	 }
	
	class PotentialAddButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			speciesDM.setSpeciesSelectionChanged(true);
			int speciesJListIndexChosen = speciesSelectionView.getSpeciesJList().getSelectedIndex();
			if(speciesDM.isSpecies1Added()){
				speciesSelectionDM.setSpeciesSelectedIndex(speciesJListIndexChosen,1);
			}
			else{
				speciesSelectionDM.setSpeciesSelectedIndex(speciesJListIndexChosen,0);
			}
			speciesSelectionView.getDescription().setText(" ");
			
			
			if(speciesJListIndexChosen == 0){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("LJ Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("LJ");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesljPotentialslist());
			}
			if(speciesJListIndexChosen == 1){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("CO2 Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("CO2");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesco2Potentialslist());
			}
			if(speciesJListIndexChosen == 4){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("methane Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("Methane");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesmethanePotentiallist());
			}
			if(speciesJListIndexChosen == 5){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("ethane Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("Ethane");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesethanePotentiallist());
			}
			if(speciesJListIndexChosen == 6){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("propane Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("Propane");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciespropanePotentiallist());
			}
			if(speciesJListIndexChosen == 7){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("alkane Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("n-Alkane");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesnalkanePotentiallist());
			}
			if(speciesJListIndexChosen == 9){
				speciesSelectionView.getDescription().append("We now create a single ");
				speciesSelectionView.getDescription().append("Water Molecule\n");
				speciesSelectionDM.setSpeciesIdentifier("Water");
				speciesSelectionDM.setPotentialList(SpeciesSelectionDM.getSpeciesh2oPotentiallist());
			}
			
			speciesSelectionView.getPotentialJList().setListData(speciesSelectionDM.getPotentialsList());
		}
		
	 }
	
	class PotentialsJListChangeListener implements ListSelectionListener
	 {
		@Override
		public void valueChanged(ListSelectionEvent e) {
			if(!speciesSelectionView.getPotentialJList().getValueIsAdjusting()){
				if(speciesDM.isSpeciesSelectionChanged()){
					if(speciesDM.isSpecies1Added()){
						speciesSelectionDM.setPotentialsSelectedIndex(
								speciesSelectionView.getPotentialJList().getSelectedIndex(), 1);
						
					}else{
						speciesSelectionDM.setPotentialsSelectedIndex(
								speciesSelectionView.getPotentialJList().getSelectedIndex(), 0);
					}
				}
			}
		}
	 }
	
	public static void main(String[] args){
		
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				SpeciesSelectionUIView s = new SpeciesSelectionUIView();
				SpeciesSelectionUIController sa = new SpeciesSelectionUIController(s,new SpeciesSelectionDM(), new SpeciesDM());
				JFrame frame = new JFrame();
				frame.add(s);
				frame.setVisible(true);
				frame.setResizable(true);
				frame.setMinimumSize(new Dimension(750,700));
				frame.setBounds(0, 0,750, 700);
				
			}
		});
		
	
	
	}
	
}
