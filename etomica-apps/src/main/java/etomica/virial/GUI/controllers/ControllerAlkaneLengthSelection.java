/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.controllers;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


import javax.swing.SwingUtilities;
import etomica.virial.GUI.models.ModelAlkaneLengthSelection;
import etomica.virial.GUI.models.ModelSelectedSpecies;
import etomica.virial.GUI.views.ViewAlkaneLengthSelection;

public class ControllerAlkaneLengthSelection {
		
	private ViewAlkaneLengthSelection[] alkaneLengthSelectionView;

	private ModelAlkaneLengthSelection modelAlkaneLengthSelection;
	private ModelSelectedSpecies modelSelectedSpecies;
	
	
	

	public ControllerAlkaneLengthSelection(ModelSelectedSpecies modelSelectedSpecies){	
		modelAlkaneLengthSelection = new ModelAlkaneLengthSelection();
		this.modelSelectedSpecies = modelSelectedSpecies;
		alkaneLengthSelectionView = new ViewAlkaneLengthSelection[10];
		
	}
	
	public ViewAlkaneLengthSelection getAlkaneLengthSelectionView(int index) {
		return alkaneLengthSelectionView[index];
	}


	public void setAlkaneLengthSelectionView(
			ViewAlkaneLengthSelection alkaneLengthSelectionView1, int index) {
		this.alkaneLengthSelectionView[index] = alkaneLengthSelectionView1;
		this.alkaneLengthSelectionView[index].openAlkaneView();
	}


	public ModelAlkaneLengthSelection getAlkaneLengthSelectionDM() {
		return modelAlkaneLengthSelection;
	}


	public void setAlkaneLengthSelectionDM(
			ModelAlkaneLengthSelection modelAlkaneLengthSelection) {
		this.modelAlkaneLengthSelection = modelAlkaneLengthSelection;
	}


	public ModelSelectedSpecies getSpeciesDM() {
		return modelSelectedSpecies;
	}


	public void setSpeciesDM(ModelSelectedSpecies modelSelectedSpecies) {
		this.modelSelectedSpecies = modelSelectedSpecies;
	}

	public void addAlkaneLengthViewListener(int index){
		getAlkaneLengthSelectionView(index).addSaveValuesButtonListener(new SaveValuesButtonListener());
		getAlkaneLengthSelectionView(index).addCloseWindowButtonListener(new CloseWindowButtonListener());
	}
	
	class SaveValuesButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			for(int i=0;i<10;i++){
				if(alkaneLengthSelectionView[i] != null){
					if(e.getSource().equals(alkaneLengthSelectionView[i].getSaveValues())){
						int Alkane1Length = Integer.parseInt(alkaneLengthSelectionView[i].getNoOfSpheres().getText());
						modelAlkaneLengthSelection.setAlkaneLength(Alkane1Length,i);
						
						alkaneLengthSelectionView[i].closeAlkaneView();
						break;
					}
				}
			} 
		}
		
		
	 }
	
	
	class CloseWindowButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			for(int i=0;i<10;i++){
				if(alkaneLengthSelectionView[i] != null){
					if(e.getSource().equals(alkaneLengthSelectionView[i].getCloseWindow())){
						alkaneLengthSelectionView[i].closeAlkaneView();
						break;
					}
				}
			}
		}
	 }
	
	public static void main(String[] args){
		
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				
				ControllerAlkaneLengthSelection sa = new ControllerAlkaneLengthSelection(new ModelSelectedSpecies());
				ViewAlkaneLengthSelection sa1 = new ViewAlkaneLengthSelection("Species1 - nAlkane");
				ViewAlkaneLengthSelection sa2 = new ViewAlkaneLengthSelection("Species1 - nAlkane");
				sa.setAlkaneLengthSelectionView(sa1,0);
				sa1.openAlkaneView();
				sa.addAlkaneLengthViewListener(0);
				//sa.setAlkaneLengthSelectionView2(sa2);
				//sa.addAlkaneLengthView2Listener();
				sa2.openAlkaneView();
				
			}
		});
		
	
	
	}
}
