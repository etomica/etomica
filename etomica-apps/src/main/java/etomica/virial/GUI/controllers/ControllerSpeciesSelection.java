/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.controllers;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import javax.swing.SwingUtilities;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableModel;

import etomica.virial.GUI.components.ICollectionPotential;

import etomica.virial.GUI.components.ExceptionDoNotRemoveSpecies;
import etomica.virial.GUI.components.ExceptionSameSpeciesSamePotential;
import etomica.virial.GUI.components.ExceptionSingleSiteAtomicPotential;
import etomica.virial.GUI.components.BuilderCollectionPotential;
import etomica.virial.GUI.components.SimulationRunner;




import etomica.virial.GUI.models.ModelSimulationEnvironment;
import etomica.virial.GUI.models.IMolecularModel_SpeciesFactory;
import etomica.virial.GUI.models.EnumPotentialParamDescription;
import etomica.virial.GUI.models.ModelSimulationConstructor;
import etomica.virial.GUI.models.ModelSelectedSpecies;
import etomica.virial.GUI.models.ModelSpeciesSelection;
import etomica.virial.GUI.models.ModelTemperatureAndSteps;
import etomica.virial.GUI.views.ViewAlertMsgBoxSpeciesRemoval;
import etomica.virial.GUI.views.ViewAlkaneLengthSelection;
import etomica.virial.GUI.views.ViewConsoleOutputBox;
import etomica.virial.GUI.views.ViewSpeciesSelection;


public class ControllerSpeciesSelection {
		
	private ViewSpeciesSelection speciesSelectionView;
	private ViewConsoleOutputBox consoleOutputBoxView;
	
	//DataModels required
	private ModelSpeciesSelection modelSpeciesSelection;
	
	private ModelSelectedSpecies modelSelectedSpecies;
	
	private ModelSimulationEnvironment simEnvDM;
	private ModelTemperatureAndSteps modelTemperatureAndSteps;
	private ModelSimulationConstructor simConstructionDM;
	
	
	//Other Controllers:
	private ControllerAlkaneLengthSelection alkaneLengthController;
	private ControllerAlertMsgBox alertMsgBoxController;
	private ControllerSimEnvironmentSelection simEnvController;
	
	
	
	ControllerSpeciesSelection(ViewSpeciesSelection speciesSelectionV, 
										ModelSpeciesSelection  modelSpeciesSelection, ModelSelectedSpecies modelSelectedSpecies ){
		this.speciesSelectionView = speciesSelectionV;
		this.consoleOutputBoxView = new ViewConsoleOutputBox();
		this.modelSpeciesSelection = modelSpeciesSelection;
		this.modelSelectedSpecies = modelSelectedSpecies;
		this.modelTemperatureAndSteps = ModelTemperatureAndSteps.getInstance();

		speciesSelectionView.addSpeciesJListListener(new SpeciesJListChangeListener());
		speciesSelectionView.addPotentialButtonListener(new PotentialAddButtonListener());
		speciesSelectionView.addPotentialsJListListener(new PotentialsJListChangeListener());
		speciesSelectionView.addEditVariablesButtonListener(new EditVariablesAddButtonListener());
		speciesSelectionView.addAlterSimEnvButtonListener(new AlterSimEnvButtonListener());
		speciesSelectionView.addRunButtonListener(new RunButtonListener());
		speciesSelectionView.setCountOfSpeciesSLMListener(new ChangeListener[10]);
		for (int i=0;i<10;i++){
			speciesSelectionView.addCountOfSpeciesSLMListener(new CountOfSpeciesJSpinnerListener(),i);
		}

		speciesSelectionView.addResetButtonListener(new ResetButtonListener());
		
	}
	
	class SpeciesJListChangeListener implements ListSelectionListener
	 {

		@Override
		public void valueChanged(ListSelectionEvent e) {
			// TODO Auto-generated method stub
			if(!speciesSelectionView.getSpeciesJList().getValueIsAdjusting()){
				
				if(!modelSpeciesSelection.isSpeciesAlreadyAdded(0)){
					modelSpeciesSelection.setMoleculeListDisplayString(new ArrayList<String>());
					speciesSelectionView.getPotentialParamDisplayTM().removeData();
					if(alkaneLengthController != null){
						alkaneLengthController = null;
					}
					
				}
				
				if(modelSpeciesSelection.getMoleculeListDisplayString().size() > 1 
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(0) 
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(1)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(2)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(3)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(4)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(5)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(6)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(7)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(8)
						&& !modelSpeciesSelection.isSpeciesAlreadyAdded(9)){
					modelSpeciesSelection.getMoleculeListDisplayString().clear();
					speciesSelectionView.getPotentialParamDisplayTM().removeData();
					speciesSelectionView.getMoleculeListDisplayTM().removeData();
				}
				speciesSelectionView.getPotentialJList().setEnabled(true);
				modelSpeciesSelection.setSpeciesSelectionChanged(true);
				
				
				
				int speciesJListIndexChosen = speciesSelectionView.getSpeciesJList().getSelectedIndex();
				//speciesSelectionView.getPotentialJList().setListData(ModelSpeciesSelection.getIntialpotentialdisplaylist());
				
				int nthSpeciesSelected = checkIndexOfSpeciesSelected();
				modelSpeciesSelection.setSpeciesIdentifier(null,nthSpeciesSelected);
				modelSpeciesSelection.setSpeciesSelectedIndex(speciesJListIndexChosen,nthSpeciesSelected);
				speciesSelectionView.getDescription().setText(" ");
				
				if(speciesJListIndexChosen == 0){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("LJ Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("LJ",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesljPotentialslist());
				}
				if(speciesJListIndexChosen == 1){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("CO2 Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("CO2",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesco2Potentialslist());
				}
				if(speciesJListIndexChosen == 2){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("methane Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("Methane",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesmethanePotentiallist());
				}
				if(speciesJListIndexChosen == 3){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("ethane Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("Ethane",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesethanePotentiallist());
				}
				if(speciesJListIndexChosen == 4){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("propane Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("Propane",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciespropanePotentiallist());
				}
				if(speciesJListIndexChosen == 5){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("alkane Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("n-Alkane",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesnalkanePotentiallist());
				}
				if(speciesJListIndexChosen == 7){
					speciesSelectionView.getDescription().append("We now create a single ");
					speciesSelectionView.getDescription().append("Water Molecule\n");
					modelSpeciesSelection.setSpeciesIdentifier("H2O",nthSpeciesSelected);
					modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesh2oPotentiallist());
				}
				
				speciesSelectionView.getPotentialJList().setListData(modelSpeciesSelection.getPotentialsList());
			}
		}
	 
	 }
	
	class PotentialAddButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			
				
			if(modelSpeciesSelection.getNthSpeciesSelectedIndex() == checkIndexOfSpeciesSelected() && modelSpeciesSelection.getNthSpeciesSelectedIndex() > 0 ||
					modelSpeciesSelection.getNthSpeciesSelectedIndex() == checkIndexOfSpeciesSelected() && modelSpeciesSelection.getNthSpeciesSelectedIndex() == 0){
				System.out.println("different molecule chosen!!");
				modelSelectedSpecies.setSpeciesDataModel(modelSpeciesSelection.getSelectedSpecies(),modelSpeciesSelection.getNthSpeciesSelectedIndex());
				modelSelectedSpecies.setNthSpeciesAdded(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
				modelSelectedSpecies.setSpecies(modelSelectedSpecies.getSpeciesDataModel(
													modelSpeciesSelection.getNthSpeciesSelectedIndex()).createSpecies(),modelSpeciesSelection.getNthSpeciesSelectedIndex());
				modelSpeciesSelection.setSpeciesAlreadyAdded(true,modelSpeciesSelection.getNthSpeciesSelectedIndex());
		
				if(modelSpeciesSelection.getSpeciesIdentifier(modelSpeciesSelection.getNthSpeciesSelectedIndex())=="n-Alkane"){
					modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).setParameter("NUMBER", 
					Integer.toString(alkaneLengthController.getAlkaneLengthSelectionDM().getAlkaneLength(modelSpeciesSelection.getNthSpeciesSelectedIndex())));
				}
			
				modelSpeciesSelection.getMoleculeListDisplayString().add(modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).
																																			getMoleculeDisplayName());
				speciesSelectionView.getJSpinnerCountOfSpecies(modelSpeciesSelection.getNthSpeciesSelectedIndex()).setEnabled(true);
				speciesSelectionView.getJSpinnerCountOfSpecies(modelSpeciesSelection.getNthSpeciesSelectedIndex()).setValue(
																												modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
			
				speciesSelectionView.addMoleculeListDisplayLSMListener(new MoleculeListDisplayLSMListener());
				speciesSelectionView.displayMoleculeList();
		
			
				
			}else{
				
				System.out.println("same molecule chosen!!");
				System.out.println("SpeciesSelectionIndex " + modelSpeciesSelection.getNthSpeciesSelectedIndex()+"\n");
				
				try {
					if(modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).getClass().getName().contains("Alkane")){
						modelSelectedSpecies.setSpeciesDataModel(modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).getClass().getConstructor(Integer.TYPE).newInstance(modelSpeciesSelection.getConstructorObject()),modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
					}else{
						modelSelectedSpecies.setSpeciesDataModel(modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).getClass().getDeclaredConstructor().newInstance(),modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
					}
				} catch (InstantiationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (SecurityException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InvocationTargetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (NoSuchMethodException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				modelSelectedSpecies.setNthSpeciesAdded(modelSpeciesSelection.getNthSpeciesSelectedIndex()+2);
				System.out.println("selectedSpeciesAddedIndex " + modelSelectedSpecies.getNthSpeciesAdded()+"\n");
				modelSelectedSpecies.setSpecies(modelSelectedSpecies.getSpeciesDataModel(
													modelSpeciesSelection.getNthSpeciesSelectedIndex()+1).createSpecies(),modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
				
		
				if(modelSpeciesSelection.getSpeciesIdentifier(modelSpeciesSelection.getNthSpeciesSelectedIndex())=="n-Alkane"){
					 alkaneLengthController.setAlkaneLengthSelectionView(new ViewAlkaneLengthSelection("Species"+ (Integer.toString(modelSpeciesSelection.getNthSpeciesSelectedIndex() + 2)) + 
								"- Alkane Chain Length"),modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
					alkaneLengthController.addAlkaneLengthViewListener(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1);
					
					modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1).setParameter("NUMBER", 
					Integer.toString(alkaneLengthController.getAlkaneLengthSelectionDM().getAlkaneLength(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1)));
					
				}
			
				modelSpeciesSelection.getMoleculeListDisplayString().add(modelSelectedSpecies.getSpeciesDataModel(modelSpeciesSelection.getNthSpeciesSelectedIndex()).
																																			getMoleculeDisplayName());
				speciesSelectionView.getJSpinnerCountOfSpecies(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1).setEnabled(true);
				speciesSelectionView.getJSpinnerCountOfSpecies(modelSpeciesSelection.getNthSpeciesSelectedIndex()+1).setValue(
																												modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
			
				speciesSelectionView.addMoleculeListDisplayLSMListener(new MoleculeListDisplayLSMListener());
				speciesSelectionView.displayMoleculeList();
		
			
				checkIndexOfSpeciesSelected();
				modelSpeciesSelection.setNthSpeciesSelectedIndex(checkIndexOfSpeciesSelected());
				System.out.println("SpeciesSelectionIndex " + modelSpeciesSelection.getNthSpeciesSelectedIndex()+"\n");
				modelSpeciesSelection.setPotentialsSelectedIndex(
						speciesSelectionView.getPotentialJList().getSelectedIndex(), modelSpeciesSelection.getNthSpeciesSelectedIndex());
				modelSpeciesSelection.setSpeciesSelectedIndex(speciesSelectionView.getSpeciesJList().getSelectedIndex(),modelSpeciesSelection.getNthSpeciesSelectedIndex());
				modelSpeciesSelection.setSpeciesIdentifier(modelSpeciesSelection.getSpeciesIdentifier(modelSpeciesSelection.getNthSpeciesSelectedIndex() - 1), modelSpeciesSelection.getNthSpeciesSelectedIndex());
				modelSpeciesSelection.setSpeciesAlreadyAdded(true,modelSpeciesSelection.getNthSpeciesSelectedIndex());
			}
			
			int totalCountMolecules = 0;
			for(int i=0;i<modelSelectedSpecies.getNthSpeciesAdded()+1;i++){
				totalCountMolecules = totalCountMolecules + Integer.parseInt((String) speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue());
			}
			speciesSelectionView.getDescription().append(Integer.toString(totalCountMolecules) + " Molecules now in sim box\n");	
			/*if(modelSpeciesSelection.getSpeciesIdentifier(modelSpeciesSelection.getNthSpeciesSelectedIndex())=="LJ"){
				speciesSelectionView.getEditVariables().setEnabled(true);
			}*/
			
		
		}
		
	 }
	
	class PotentialsJListChangeListener implements ListSelectionListener
	 {
		@SuppressWarnings("unchecked")
		@Override
		public void valueChanged(ListSelectionEvent e) {
			if(!speciesSelectionView.getPotentialJList().getValueIsAdjusting()){
				if(modelSpeciesSelection.isSpeciesSelectionChanged()){
					
					int currentPotentialIndex = 0;
					
					currentPotentialIndex = speciesSelectionView.getPotentialJList().getSelectedIndex();
					
					if(currentPotentialIndex != -1){
						
						
						//This has use for alkane contructors where you pass a Int Object
						//Object[] temporaryIntObject = new Object[1];
					
					
						//This if condition depicts if We re adding a new or a second species 
						modelSpeciesSelection.setNthSpeciesSelectedIndex(checkIndexOfSpeciesSelected());
						int nthSpeciesSelected = modelSpeciesSelection.getNthSpeciesSelectedIndex();
						modelSpeciesSelection.setPotentialsSelectedIndex(
							speciesSelectionView.getPotentialJList().getSelectedIndex(), nthSpeciesSelected);
						
					
						speciesSelectionView.getAddPotential().setEnabled(true);
						//speciesSelectionView.getPotentialParamDispayTable().setEnabled(true);
						speciesSelectionView.getMoleculeListDisplayTable().setEnabled(true);
					
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="LJ"){
							speciesSelectionView.getEditVariables().setEnabled(true);
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single LJ Molecule\n");
								speciesSelectionView.getDescription().append("With a Spherical-2-Body Potential\n");
							}else if(currentPotentialIndex == 1){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single LJ Molecule\n");
								speciesSelectionView.getDescription().append("With a Spherical-2-Body-With-Quad Potential\n");
							
							}else if(currentPotentialIndex == 2){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single LJ Molecule\n");
								speciesSelectionView.getDescription().append("With a 2-Centred-With-Quad Potential\n");
							}
						
						}
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="CO2"){
							speciesSelectionView.getEditVariables().setEnabled(false);
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single CO2 Molecule\n");
								speciesSelectionView.getDescription().append("With a 2-Centred-With-Quad Potential\n");
							}
					
						if(currentPotentialIndex == 1){
							speciesSelectionView.getDescription().setText(" ");
							speciesSelectionView.getDescription().append("We now create a single CO2 Molecule\n");
							speciesSelectionView.getDescription().append("Modeled as EPM2\n");
							
							}
						if(currentPotentialIndex == 2){
							speciesSelectionView.getDescription().setText(" ");
							speciesSelectionView.getDescription().append("We now create a single CO2 Molecule\n");
							speciesSelectionView.getDescription().append("Modeled as TRAPPE\n");
							}
						}
					
					
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="Methane"){
							speciesSelectionView.getEditVariables().setEnabled(false);
							//temporaryIntObject[0] = new Integer(1);
							modelSpeciesSelection.setConstructorObject(new Integer(1));
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single Methane Molecule\n");
								speciesSelectionView.getDescription().append("Modelled as TRAPPE\n");
							}
						}
					
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="Ethane"){
							speciesSelectionView.getEditVariables().setEnabled(false);
							//temporaryIntObject[0] = new Integer(2);
							modelSpeciesSelection.setConstructorObject(new Integer(2));
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single Ethane Molecule\n");
								speciesSelectionView.getDescription().append("Modelled as TRAPPE\n");
							}
						
						}
					
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="Propane"){
							speciesSelectionView.getEditVariables().setEnabled(false);
							//temporaryIntObject[0] = new Integer(3);
							modelSpeciesSelection.setConstructorObject(new Integer(3));
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single Propane Molecule\n");
								speciesSelectionView.getDescription().append("Modelled as TRAPPE\n");
							}
						
						}
					
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="n-Alkane"){
							speciesSelectionView.getEditVariables().setEnabled(false);
						
							//temporaryIntObject[0] = new Integer(4);
							modelSpeciesSelection.setConstructorObject(new Integer(4));
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single higher alkane Molecule\n");
								speciesSelectionView.getDescription().append("Modelled as TRAPPE\n");
							}
							 if(alkaneLengthController == null){
								 alkaneLengthController = new ControllerAlkaneLengthSelection(modelSelectedSpecies);
								 alkaneLengthController.setAlkaneLengthSelectionView(new ViewAlkaneLengthSelection("Species"+ (Integer.toString(modelSpeciesSelection.getNthSpeciesSelectedIndex() + 1)) + 
										"- Alkane Chain Length"),modelSpeciesSelection.getNthSpeciesSelectedIndex());
								 alkaneLengthController.addAlkaneLengthViewListener(modelSpeciesSelection.getNthSpeciesSelectedIndex());
							}
							 else{
								 alkaneLengthController.setAlkaneLengthSelectionView(new ViewAlkaneLengthSelection("Species"+ (Integer.toString(modelSpeciesSelection.getNthSpeciesSelectedIndex() + 1)) + 
											"- Alkane Chain Length"),modelSpeciesSelection.getNthSpeciesSelectedIndex());
								 alkaneLengthController.addAlkaneLengthViewListener(modelSpeciesSelection.getNthSpeciesSelectedIndex());
							 }
						}
						if(modelSpeciesSelection.getSpeciesIdentifier(nthSpeciesSelected)=="H2O"){
							speciesSelectionView.getEditVariables().setEnabled(false);
							if(currentPotentialIndex == 0){
								speciesSelectionView.getDescription().setText(" ");
								speciesSelectionView.getDescription().append("We now create a single Water Molecule\n");
								speciesSelectionView.getDescription().append("Modelled as SPC\n");
							}
						}
					
					
						if(!speciesSelectionView.getReset().isEnabled()){
							speciesSelectionView.getReset().setEnabled(true);
						}
						if(!speciesSelectionView.getRun().isEnabled()){
							speciesSelectionView.getRun().setEnabled(true);
						}
						if(!speciesSelectionView.getAlterSimEnv().isEnabled()){
							speciesSelectionView.getAlterSimEnv().setEnabled(true);
						}
					
						try{
												

							@SuppressWarnings("rawtypes")
							Constructor[] speciesFactoryConstructor = modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].getConstructors();
							if(speciesFactoryConstructor.length > 1){
								for(int i = 0; i < speciesFactoryConstructor.length;i++){
									if(speciesFactoryConstructor[i].getParameterTypes().length > 0){
									
										if(nthSpeciesSelected != 0){
											IMolecularModel_SpeciesFactory tempSpecies3ClassObject = (IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].getConstructor(Integer.TYPE).newInstance(modelSpeciesSelection.getConstructorObject());
										
											for(int m=0;m<nthSpeciesSelected;m++){
												if(checkIfSameSpeciesWithSimilarPotentials(tempSpecies3ClassObject,modelSelectedSpecies.getSpeciesDataModel(m),
																						modelSpeciesSelection.getSpeciesSelectedIndex(nthSpeciesSelected),
																						modelSpeciesSelection.getSpeciesSelectedIndex(m))){
													setAlertMsgBoxController(new ControllerAlertMsgBox("You ve chosen the same potential. Please choose another same",modelSpeciesSelection));
													throw new ExceptionSameSpeciesSamePotential();
												}
											}
										
											modelSpeciesSelection.setSelectedSpecies((IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].
													getConstructor(Integer.TYPE).newInstance(modelSpeciesSelection.getConstructorObject()));
											
									
										}
										else{
										
											IMolecularModel_SpeciesFactory tempSpeciesClassObject = (IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].
																										getConstructor(Integer.TYPE).newInstance(modelSpeciesSelection.getConstructorObject());
										
											for(int m=0;m<nthSpeciesSelected;m++){
												if(checkIfSameSpeciesWithSimilarPotentials(tempSpeciesClassObject,modelSelectedSpecies.getSpeciesDataModel(m),
																						modelSpeciesSelection.getSpeciesSelectedIndex(nthSpeciesSelected),
																						modelSpeciesSelection.getSpeciesSelectedIndex(m))){
													setAlertMsgBoxController(new ControllerAlertMsgBox("You ve chosen the same potential. Please choose another same",modelSpeciesSelection));
													throw new ExceptionSameSpeciesSamePotential();
												}
											}
										
											modelSpeciesSelection.setSelectedSpecies((IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].
													getConstructor(Integer.TYPE).newInstance(modelSpeciesSelection.getConstructorObject()));
										}
									
									}
									
								}
							}else{
							
								IMolecularModel_SpeciesFactory tempSpeciesClassObject = (IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].getConstructor().newInstance(new Object[0]);
							
							
								for(int m=0;m<nthSpeciesSelected;m++){
									if(checkIfSameSpeciesWithSimilarPotentials(tempSpeciesClassObject,modelSelectedSpecies.getSpeciesDataModel(m),
																			modelSpeciesSelection.getSpeciesSelectedIndex(nthSpeciesSelected),
																			modelSpeciesSelection.getSpeciesSelectedIndex(m))){
										setAlertMsgBoxController(new ControllerAlertMsgBox("You ve chosen the same potential. Please choose another same",modelSpeciesSelection));
										throw new ExceptionSameSpeciesSamePotential();
									}
								}

								modelSpeciesSelection.setSelectedSpecies((IMolecularModel_SpeciesFactory)modelSpeciesSelection.getPotentialsClassList()[currentPotentialIndex].getConstructor().newInstance(new Object[0]));
							}
						
							modelSpeciesSelection.setPotentialParamString(modelSpeciesSelection.getSelectedSpecies().getParametersArray());
							modelSpeciesSelection.setPotentialParamAndValuesString(modelSpeciesSelection.getSelectedSpecies().getParamAndValues());
							modelSpeciesSelection.setSitesOnSpecies(modelSpeciesSelection.getSelectedSpecies().getPotentialSites());
						
							if(!speciesSelectionView.getPotentialParamDisplayLSM().isSelectionEmpty()){
								speciesSelectionView.getPotentialParamDisplayLSM().clearSelection();
							}
							speciesSelectionView.getPotentialParamDisplayTM().removeData();
							for(int i = 0;i <modelSpeciesSelection.getPotentialParamAndValuesString().length;i++){
								for(EnumPotentialParamDescription parameters : EnumPotentialParamDescription.values()){
									if(modelSpeciesSelection.getPotentialParamAndValuesString()[i][0].toUpperCase().contains(parameters.toString())){
										speciesSelectionView.getPotentialParamDisplayTM().setValueAt(parameters.unit() != null ? modelSpeciesSelection.getPotentialParamAndValuesString()[i][0]+" "+parameters.unit() : modelSpeciesSelection.getPotentialParamAndValuesString()[i][0],i,0);
										speciesSelectionView.getPotentialParamDisplayTM().setValueAt(modelSpeciesSelection.getPotentialParamAndValuesString()[i][1],i,1);
									}
								}
							}
							
						}
						catch (ExceptionSameSpeciesSamePotential E){
							resetViewAfterSSSPException();
						}
						catch(Exception refException){
							refException.printStackTrace();
						}
					
					
					}else{
						speciesSelectionView.getAddPotential().setEnabled(false);
					}
					
				}
			}
		}
		
		
	 }
	
	
	
	
	class EditVariablesAddButtonListener implements ActionListener 
	 {
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			speciesSelectionView.getPotentialParamDispayTable().setEnabled(true);
			speciesSelectionView.addPotParamDisplayTMListener(new PotParamDisplayTMListener());
			//speciesSelectionView.getMoleculeListDisplayTable().setEnabled(true);	
			//speciesSelectionView.addMoleculeListDisplayTMListener(new MoleculeListDisplayTMListener());
		}

	 }
	
	class PotParamDisplayTMListener implements TableModelListener 
	 {

		@Override
		public void tableChanged(TableModelEvent e) {
			Object dataValue= null;
			Object dataStringWithUnit= null;
			StringTokenizer dataString = null;
			int row = e.getFirstRow();
	        int column = e.getColumn();
	        int speciesListIndex = 0;
	        TableModel model = (TableModel)e.getSource();
	        try{
	        	if(column != -1){
	        		dataValue = model.getValueAt(row, column);
	        		dataStringWithUnit = model.getValueAt(row, column - 1);
	        		dataString = new StringTokenizer((String)dataStringWithUnit,"(");
	        		if(speciesSelectionView.getMoleculeListDisplayTable().getSelectedRows().length != 0){
		        		int[] selectedRows = speciesSelectionView.getMoleculeListDisplayTable().getSelectedRows();
		        		speciesListIndex = selectedRows[0];
		        		System.out.println(speciesListIndex);
		        	}
		        	if(modelSpeciesSelection.getSpeciesIdentifier(speciesListIndex) == "LJ"){
		        			speciesSelectionView.getDescription().append("You have now modified the default value!");
		        			modelSelectedSpecies.getSpeciesDataModel(speciesListIndex).setParameter(dataString.nextToken().trim(), (String)dataValue);
		        			speciesSelectionView.getPotentialParamDispayTable().clearSelection();
		        			
		        			speciesSelectionView.getPotentialParamDispayTable().setEnabled(false);
		        			speciesSelectionView.getEditVariables().setEnabled(false);
		        			
		        			if(speciesSelectionView.getMoleculeListDisplayTable().getSelectedRows().length != 0){
		        				speciesSelectionView.getMoleculeListDisplayTable().clearSelection();
		        			}
		        			
		        	}
	        	}
	        	
	        }
	        catch(NullPointerException n){
	        		
	        }
		}
	 }
	
	class MoleculeListDisplayTMListener implements TableModelListener 
	 {

		@Override
		public void tableChanged(TableModelEvent e) {
			// TODO Auto-generated method stub
			
	        try{
	        	System.out.println("HI");
	        	
	        }
	        catch(NullPointerException n){
	        		
	        }
		}
	
	 }
	
	
	class PotParamDisplayLSMListener implements ListSelectionListener 
	 {

		@Override
		public void valueChanged(ListSelectionEvent e) {
			// TODO Auto-generated method stub
			if( !e.getValueIsAdjusting()){
				
				String selectedParameter = null;
				int RowIndex = 0;
				if(speciesSelectionView.getPotentialParamDispayTable().getSelectedRows().length != 0){
					int[] selectedRows = speciesSelectionView.getPotentialParamDispayTable().getSelectedRows();
					RowIndex = selectedRows[0];
				}
				try{
					selectedParameter = (String)speciesSelectionView.getPotentialParamDispayTable().getValueAt(RowIndex,0);
					for(EnumPotentialParamDescription parameters : EnumPotentialParamDescription.values()){
						if(selectedParameter.toUpperCase().contains(parameters.toString())){
							if(speciesSelectionView.getEditVariables().isEnabled()){
								speciesSelectionView.getDescription().append("You have now chosen to edit the default value of " + parameters.toString() + "\n");
							}
							speciesSelectionView.getDescription().append("This chosen parameter is " + parameters.description()+"\n");
						}
					}
				}
				catch(NullPointerException n){}
			}
		}
	 }
	
	
	class MoleculeListDisplayLSMListener implements ListSelectionListener 
	 {
		@Override
		public void valueChanged(ListSelectionEvent e) {
			// TODO Auto-generated method stub
			
			if( !e.getValueIsAdjusting()){
				System.out.println("there");
				speciesSelectionView.getPotentialParamDispayTable().clearSelection();
				speciesSelectionView.removePotParamDisplayTMListener();
				String selectedParameter = null;
				String[][] ParamAndValue = null;
				int RowIndex = 0;
				if(speciesSelectionView.getMoleculeListDisplayTable().getSelectedRows().length != 0){
					int[] selectedRows = speciesSelectionView.getMoleculeListDisplayTable().getSelectedRows();
					RowIndex = selectedRows[0];
					try{
						selectedParameter = (String)speciesSelectionView.getMoleculeListDisplayTable().getValueAt(RowIndex,1);
						if(selectedParameter == modelSelectedSpecies.getSpeciesDataModel(RowIndex).getMoleculeDisplayName()){
							ParamAndValue = modelSelectedSpecies.getSpeciesDataModel(RowIndex).getParamAndValues();
							if(modelSpeciesSelection.getSpeciesIdentifier(RowIndex)!="LJ"){
								speciesSelectionView.getEditVariables().setEnabled(false);
							}else{
								speciesSelectionView.getEditVariables().setEnabled(true);
							}
						}
						speciesSelectionView.getPotentialParamDisplayTM().removeData();
						for(int i = 0;i <ParamAndValue.length;i++){
							for(EnumPotentialParamDescription parameters : EnumPotentialParamDescription.values()){
								if(ParamAndValue[i][0].toUpperCase().contains(parameters.toString())){
									speciesSelectionView.getPotentialParamDisplayTM().setValueAt(parameters.unit() != null ? ParamAndValue[i][0]+" "+parameters.unit() : ParamAndValue[i][0],i,0);
									speciesSelectionView.getPotentialParamDisplayTM().setValueAt(ParamAndValue[i][1],i,1);
								}
							}
						}
					}
					catch(NullPointerException n){}
				}
				else{
					speciesSelectionView.getEditVariables().setEnabled(false);
				}
				
			}
		}
	 }
	
	
	class AlterSimEnvButtonListener implements ActionListener 
	 {
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			if(simEnvDM == null){
				if(modelSelectedSpecies.getSpeciesDataModel(0) == null){
					simEnvDM = new ModelSimulationEnvironment(modelTemperatureAndSteps.getTemperature(),modelTemperatureAndSteps.getNoOfSteps());
				}
				else{
					simEnvDM = new ModelSimulationEnvironment(modelTemperatureAndSteps.getTemperature(),modelTemperatureAndSteps.getNoOfSteps(),
																											modelSelectedSpecies);
				}
			}else{
				
				int nthSpeciesSelected = checkIndexOfSpeciesSelected();
				for (int j=0;j<modelSelectedSpecies.getNthSpeciesAdded();j++){
					if(modelSelectedSpecies.getSpeciesDataModel(j) != null){
						simEnvDM.setSigmaHSRef(simEnvDM.getSigmaHSRefSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
						simEnvDM.setLengthAlkaneChain(simEnvDM.getAlkaneSpheresSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
						//simEnvDM.setCountSpecies(modelSelectedSpecies.getSpeciesMoleculeCount(j),j);
					}
				}
			}
			
			simEnvController = new ControllerSimEnvironmentSelection(simEnvDM);
		}

	 }
	
	class RunButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			simConstructionDM = new ModelSimulationConstructor();
			//simConstructionDM.setFactory(new SimulationRunner(modelSelectedSpecies));
			if(modelTemperatureAndSteps == null){
				modelTemperatureAndSteps = ModelTemperatureAndSteps.getInstance();
			}
			if(simEnvDM == null){
				simEnvDM = new ModelSimulationEnvironment(modelTemperatureAndSteps.getTemperature(),modelTemperatureAndSteps.getNoOfSteps(),
																											modelSelectedSpecies);
			}
			
			for (int j=0;j<modelSelectedSpecies.getNthSpeciesAdded();j++){
				if(modelSelectedSpecies.getSpeciesDataModel(j) != null){
					simEnvDM.setSigmaHSRef(simEnvDM.getSigmaHSRefSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
					simEnvDM.setLengthAlkaneChain(simEnvDM.getAlkaneSpheresSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
					//simEnvDM.setCountSpecies(modelSelectedSpecies.getSpeciesMoleculeCount(j),j);
				}
			}

			simEnvDM.setnPoints(modelSelectedSpecies.getSpeciesMoleculeCount());
			
			
			simEnvDM.setSystemSigmaHSRef(simEnvDM.calculateSystemSigmaHSRef(modelSelectedSpecies.getSpeciesDataModel(),modelSelectedSpecies.getSpeciesMoleculeCount()));
			
			
			simEnvDM.calculateSystemSigmaHSRef(modelSelectedSpecies);
			simConstructionDM.setPotentialsCollectionBuilder(new BuilderCollectionPotential());
			try{
				ArrayList<ICollectionPotential> tempPotentialCollections = 
						simConstructionDM.getPotentialsCollectionBuilder().checkIfCompatible(modelSelectedSpecies,simEnvDM );
				if(tempPotentialCollections == null){
					if(simEnvDM.isStopSimulation() == false){
						throw new ExceptionSingleSiteAtomicPotential();
					}else{
						setAlertMsgBoxController(new ControllerAlertMsgBox("Oops,We cannot run this simulation",modelSpeciesSelection));
						resetSelectionViewANew();
					}
					
				}else{
					simConstructionDM.setArrayListPotentialCollection(tempPotentialCollections);
					SimulationRunner simRun  = new SimulationRunner(modelSelectedSpecies);
					simConstructionDM.setFactory(simRun);
					simConstructionDM.getFactory().runSimulation(simEnvDM, tempPotentialCollections);
				}
			}
			catch(ExceptionSingleSiteAtomicPotential a){
				//simConstructionDM.setFactory(new SimulationRunner(modelSelectedSpecies));
				//simConstructionDM.getFactory().runSimulation(simEnvDM);
			}catch(Exception exception){
				setAlertMsgBoxController(new ControllerAlertMsgBox("Oops, We've encountered some problem!",modelSpeciesSelection));
				exception.printStackTrace();
			}
			
		}
	 }
	
	
	class ResetButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent e) {
			resetSelectionViewANew();
			
			for(int i=0;i<10;i++){
				speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[0]);
				speciesSelectionView.getJSpinnerCountOfSpecies(i).setEnabled(false);
			}
			 
		}
		
	 }	
	

	class CountOfSpeciesJSpinnerListener implements ChangeListener 
	 {

		@SuppressWarnings("null")
		@Override
		public void stateChanged(ChangeEvent event) {
			// TODO Auto-generated method stub
			
			int totalCountMolecules = 0;
			for(int i=0;i<modelSelectedSpecies.getNthSpeciesAdded();i++){
				totalCountMolecules = totalCountMolecules + Integer.parseInt((String) speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue());
			}
			speciesSelectionView.getDescription().append(Integer.toString(totalCountMolecules) + " Molecules now in sim box\n");
			if(totalCountMolecules > 9){
				setAlertMsgBoxController(new ControllerAlertMsgBox("Please reduce no of molecules in any of the chosen species!!",modelSpeciesSelection));
			}
			else{
				try{
					for(int i=0;i<modelSelectedSpecies.getNthSpeciesAdded();i++){
						if(event.getSource().equals(speciesSelectionView.getCountOfSpeciesSLM(i))){
							
							if(Integer.parseInt((String)speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue()) == 0 &&
								i == 0 && i == modelSelectedSpecies.getNthSpeciesAdded() - 1){
								//last molecule removed
								
								int response = JOptionPane.showConfirmDialog(null, "Do you want to remove species "+ Integer.toString(i+1)+"?", "Confirm",
								        JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);

								if(response == JOptionPane.NO_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}else if(response == JOptionPane.YES_OPTION){
									resetSelectionViewANew();
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setEnabled(false);
									
									break;
								}else if(response == JOptionPane.CLOSED_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}
								
								
							}else if(Integer.parseInt((String)speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue()) == 0 &&
										i < (modelSelectedSpecies.getNthSpeciesAdded() - 1)){
								int response = JOptionPane.showConfirmDialog(null, "Do you want to remove species "+ Integer.toString(i+1)+"?", "Confirm",
								        JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
								
								if(response == JOptionPane.NO_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}else if(response == JOptionPane.YES_OPTION){
									modelSpeciesSelection.getMoleculeListDisplayString().remove(i);
									for(int k=i;k<modelSelectedSpecies.getNthSpeciesAdded()-1;k++){
										modelSelectedSpecies.setSpeciesDataModel(modelSelectedSpecies.getSpeciesDataModel(k+1),k);
										modelSelectedSpecies.setSpecies(modelSelectedSpecies.getSpecies(k+1),k);
										modelSelectedSpecies.setSpeciesMoleculeCount(Integer.parseInt((String) speciesSelectionView.getJSpinnerCountOfSpecies(k+1).getValue()),k);
									
										modelSpeciesSelection.setSpeciesSelectedIndex(modelSpeciesSelection.getSpeciesSelectedIndex(k+1), k);
										modelSpeciesSelection.setPotentialsSelectedIndex(modelSpeciesSelection.getPotentialsSelectedIndex()[k+1], k);
										
										speciesSelectionView.getJSpinnerCountOfSpecies(k).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[modelSelectedSpecies.getSpeciesMoleculeCount(k)]);
										
									}
									
									speciesSelectionView.removeCountOfSpeciesSLMListener(modelSelectedSpecies.getNthSpeciesAdded()-1);
									speciesSelectionView.getJSpinnerCountOfSpecies(modelSelectedSpecies.getNthSpeciesAdded()-1).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[0]);
									speciesSelectionView.getJSpinnerCountOfSpecies(modelSelectedSpecies.getNthSpeciesAdded()-1).setEnabled(false);
									speciesSelectionView.addCountOfSpeciesSLMListener(new CountOfSpeciesJSpinnerListener(),modelSelectedSpecies.getNthSpeciesAdded()-1);
									
									modelSelectedSpecies.setSpeciesDataModel(null,modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSelectedSpecies.setSpecies(null,modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSelectedSpecies.setSpeciesMoleculeCount(0,modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSpeciesSelection.setPotentialsSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSpeciesSelection.setSpeciesSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSpeciesSelection.setSpeciesAlreadyAdded(false,modelSelectedSpecies.getNthSpeciesAdded()-1);
									modelSpeciesSelection.setPotentialParamString(modelSelectedSpecies.getSpeciesDataModel(0).getParametersArray());
									modelSpeciesSelection.setPotentialParamAndValuesString(modelSelectedSpecies.getSpeciesDataModel(0).getParamAndValues());
									
									
									displaySelectedSpecies(modelSpeciesSelection.getSpeciesSelectedIndex(0),modelSpeciesSelection.getPotentialsSelectedIndex()[0]);
									speciesSelectionView.getDescription().setText("");
									modelSelectedSpecies.setNthSpeciesAdded(modelSelectedSpecies.getNthSpeciesAdded() - 1);
									System.out.println(modelSelectedSpecies.getNthSpeciesAdded());
									
									break;
								}else if(response == JOptionPane.CLOSED_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}

							}else if(Integer.parseInt((String)speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue()) == 0 &&
									i == (modelSelectedSpecies.getNthSpeciesAdded() - 1)){
								
								int response = JOptionPane.showConfirmDialog(null, "Do you want to remove species "+ Integer.toString(i+1)+"?", "Confirm",
								        JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);

								if(response == JOptionPane.NO_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}else if(response == JOptionPane.YES_OPTION){
									modelSelectedSpecies.setSpeciesDataModel(null,modelSelectedSpecies.getNthSpeciesAdded() - 1);
									modelSelectedSpecies.setSpecies(null,modelSelectedSpecies.getNthSpeciesAdded() - 1);
									modelSelectedSpecies.setSpeciesMoleculeCount(0,modelSelectedSpecies.getNthSpeciesAdded() - 1);
									
									modelSpeciesSelection.setSpeciesAlreadyAdded(false,modelSelectedSpecies.getNthSpeciesAdded() - 1);
									modelSpeciesSelection.getMoleculeListDisplayString().remove(modelSelectedSpecies.getNthSpeciesAdded() - 1);
									
									modelSpeciesSelection.setSpeciesSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded() - 1);
									modelSpeciesSelection.setPotentialsSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded() - 1);
									
									speciesSelectionView.getJSpinnerCountOfSpecies(modelSelectedSpecies.getNthSpeciesAdded() - 1).setEnabled(false);
							
									modelSpeciesSelection.setPotentialParamString(modelSelectedSpecies.getSpeciesDataModel(0).getParametersArray());
									modelSpeciesSelection.setPotentialParamAndValuesString(modelSelectedSpecies.getSpeciesDataModel(0).getParamAndValues());
									displaySelectedSpecies(modelSpeciesSelection.getSpeciesSelectedIndex(0),modelSpeciesSelection.getPotentialsSelectedIndex()[0]);
									modelSelectedSpecies.setNthSpeciesAdded(modelSelectedSpecies.getNthSpeciesAdded() - 1);
									System.out.println(modelSelectedSpecies.getNthSpeciesAdded());
									break;
								
								}else if(response == JOptionPane.CLOSED_OPTION){
									speciesSelectionView.getJSpinnerCountOfSpecies(i).setValue(modelSpeciesSelection.getCountOfMoleculesStrings()[1]);
									break;
								}

							}
							else{
								modelSelectedSpecies.setSpeciesMoleculeCount(Integer.parseInt((String) speciesSelectionView.getJSpinnerCountOfSpecies(i).getValue()),i);
								break;
							}
						}
					}
				}
				catch(Exception exception){

				}
			}
		}		
	 }
		
	public ControllerAlertMsgBox getAlertMsgBoxController() {
		return alertMsgBoxController;
	}


	public void setAlertMsgBoxController(
			ControllerAlertMsgBox alertMsgBoxController) {
		this.alertMsgBoxController = alertMsgBoxController;
	}


	public ControllerSimEnvironmentSelection getSimEnvController() {
		return simEnvController;
	}


	public void setSimEnvController(
			ControllerSimEnvironmentSelection simEnvController) {
		this.simEnvController = simEnvController;
	}
	
	public void resetSelectionViewANew(){
		modelSelectedSpecies.reset();
		modelSpeciesSelection.reset();
		simEnvDM = null;
		//speciesSelectionView.initComponents();
		//reset();
		
		speciesSelectionView.removeSpeciesJListListener();
		speciesSelectionView.getSpeciesJList().clearSelection();
		//speciesSelectionView.getSpeciesJList().setListData(ModelSpeciesSelection.getIntialdisplaylist());
		
		speciesSelectionView.getSpeciesJList().setListData(ModelSpeciesSelection.getSpeciesList());
		
		speciesSelectionView.removePotParamDisplayLSMListener();
		speciesSelectionView.getPotentialParamDisplayLSM().clearSelection();
		
		speciesSelectionView.removePotentialJListListener();
		speciesSelectionView.getPotentialJList().clearSelection();
		speciesSelectionView.getPotentialJList().setListData(ModelSpeciesSelection.getIntialpotentialdisplaylist());
		
		speciesSelectionView.getReset().setEnabled(false);
		speciesSelectionView.getRun().setEnabled(false);
		speciesSelectionView.getAlterSimEnv().setEnabled(false);
		
		speciesSelectionView.removePotParamDisplayTMListener();
		speciesSelectionView.getPotentialParamDisplayTM().removeData();
		
		speciesSelectionView.removeMoleculeListDisplayTMListener();
		speciesSelectionView.getMoleculeListDisplayTM().removeData();

		
		speciesSelectionView.addSpeciesJListListener(new SpeciesJListChangeListener());	
		speciesSelectionView.addPotentialsJListListener(new PotentialsJListChangeListener());
		
		speciesSelectionView.addPotParamDisplayLSMListener(new PotParamDisplayLSMListener());
		speciesSelectionView.addMoleculeListDisplayLSMListener(new MoleculeListDisplayLSMListener());
		speciesSelectionView.getDescription().append("Last molecule was removed!!");
		//speciesSelectionView.addPotParamDisplayTMListener(new PotParamDisplayTMListener());
		//modelSpeciesSelection.setSpeciesAlreadyAdded(false);
		
		
	}
	
	public void resetViewToDisplaySpecies1(){
		displaySelectedSpecies(modelSpeciesSelection.getSpeciesSelectedIndex(0),modelSpeciesSelection.getPotentialsSelectedIndex()[0]);		
		speciesSelectionView.getAddPotential().setEnabled(false);
		speciesSelectionView.getEditVariables().setEnabled(false);
		modelSpeciesSelection.setSpeciesAlreadyAdded(false,1);	
	}
	
	
	public void resetViewAfterSSSPException(){
		int selectedSpeciesIndex = 0;
		int selectedPotentialIndex = 0;
		
		speciesSelectionView.getJSpinnerCountOfSpecies(modelSelectedSpecies.getNthSpeciesAdded()).setEnabled(false);
		modelSpeciesSelection.setSpeciesSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded());
		modelSpeciesSelection.setPotentialsSelectedIndex(-1, modelSelectedSpecies.getNthSpeciesAdded());
		
		selectedSpeciesIndex = modelSpeciesSelection.getSpeciesSelectedIndex(modelSelectedSpecies.getNthSpeciesAdded() - 1);
		selectedPotentialIndex = modelSpeciesSelection.getPotentialsSelectedIndex()[modelSelectedSpecies.getNthSpeciesAdded() - 1];
		displaySelectedSpecies(selectedSpeciesIndex,selectedPotentialIndex);
		
	}
	
	
	
	public void displaySelectedSpecies(int speciesIndex, int potentialIndex){
		
		speciesSelectionView.removeSpeciesJListListener();
		speciesSelectionView.getSpeciesJList().clearSelection();
		speciesSelectionView.getSpeciesJList().setListData(ModelSpeciesSelection.getSpeciesList());
		speciesSelectionView.getSpeciesJList().setSelectedIndex(speciesIndex);
		
		speciesSelectionView.removePotParamDisplayLSMListener();
		speciesSelectionView.getPotentialParamDisplayLSM().clearSelection();
		
		speciesSelectionView.removePotentialJListListener();
		speciesSelectionView.getPotentialJList().clearSelection();
		
		switch(speciesIndex){
		case 0:
			//modelSpeciesSelection.setSpeciesIdentifier("LJ",nthSpeciesIndex);
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesljPotentialslist());
			
			break;
		case 1:
			//modelSpeciesSelection.setSpeciesIdentifier("CO2",nthSpeciesIndex);
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesco2Potentialslist());
			break;
		case 2:
			//modelSpeciesSelection.setSpeciesIdentifier("Methanol",nthSpeciesIndex);
			
			break;
		case 3:
			//modelSpeciesSelection.setSpeciesIdentifier("Ethanol");
			break;
		case 4:
			//modelSpeciesSelection.setSpeciesIdentifier("Methane");
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesmethanePotentiallist());
			break;
			
		case 5:
			//modelSpeciesSelection.setSpeciesIdentifier("Ethane");
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesethanePotentiallist());
			break;
		case 6:
			//modelSpeciesSelection.setSpeciesIdentifier("Propane");
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciespropanePotentiallist());
			break;
		case 7:
			//modelSpeciesSelection.setSpeciesIdentifier("Napthalene");
			break;
		case 9:
			//modelSpeciesSelection.setSpeciesIdentifier("Water");
			modelSpeciesSelection.setPotentialList(ModelSpeciesSelection.getSpeciesh2oPotentiallist());
			break;
		}
		speciesSelectionView.getPotentialJList().setListData(modelSpeciesSelection.getPotentialsList());
		speciesSelectionView.getPotentialJList().setSelectedIndex(potentialIndex);
		
		
		speciesSelectionView.removePotParamDisplayTMListener();
		speciesSelectionView.getPotentialParamDisplayTM().removeData();
		speciesSelectionView.addSpeciesJListListener(new SpeciesJListChangeListener());	

		speciesSelectionView.addPotentialsJListListener(new PotentialsJListChangeListener());
		
		speciesSelectionView.addPotParamDisplayLSMListener(new PotParamDisplayLSMListener());
		
		for(int i = 0;i <modelSpeciesSelection.getPotentialParamAndValuesString().length;i++){
			speciesSelectionView.getPotentialParamDisplayTM().setValueAt(modelSpeciesSelection.getPotentialParamAndValuesString()[i][0],i,0);
			speciesSelectionView.getPotentialParamDisplayTM().setValueAt(modelSpeciesSelection.getPotentialParamAndValuesString()[i][1],i,1);
		}
		
		speciesSelectionView.displayMoleculeList();
		//speciesSelectionView.addPotParamDisplayTMListener(new PotParamDisplayTMListener());
	}
	
	
	public boolean checkIfSameSpeciesWithSimilarPotentials(IMolecularModel_SpeciesFactory potential1, IMolecularModel_SpeciesFactory potential2, 
			int species1JListSelectedIndex,int species2JListSelectedIndex){
		if(potential1.getMoleculeDisplayName() == potential2.getMoleculeDisplayName() && alkaneLengthController != null){
			return false;
		}else if(species1JListSelectedIndex == species2JListSelectedIndex && 
				potential1.getMoleculeDisplayName() != potential2.getMoleculeDisplayName()){
			return true;
		}else{
			return false;
		}
		
	}
	
	public int checkIndexOfSpeciesSelected(){
		int index=0;
		for(int i=0;i<10;i++){
			if(!modelSpeciesSelection.isSpeciesAlreadyAdded(i)){
				index = i;
				break;
			}
		}
		return index;
	}
	public static void main(String[] args){
		
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				ModelSpeciesSelection modelSpeciesSelection = new ModelSpeciesSelection();
				ModelSelectedSpecies modelSelectedSpecies = new ModelSelectedSpecies();
				ViewSpeciesSelection s = new ViewSpeciesSelection(modelSpeciesSelection,modelSelectedSpecies);
				ControllerSpeciesSelection sa = new ControllerSpeciesSelection(s,modelSpeciesSelection,modelSelectedSpecies);
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
