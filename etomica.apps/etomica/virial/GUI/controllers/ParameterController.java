package etomica.virial.GUI.controllers;

import java.awt.CardLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JComboBox;
import javax.swing.JOptionPane;

import etomica.virial.GUI.containers.SingleSpeciesParam;
import etomica.virial.GUI.models.LJDefaultParametersModel;
import etomica.virial.GUI.models.SingleSpeciesModel;
import etomica.virial.GUI.models.SpeciesModel;
import etomica.virial.GUI.models.SuperModel;
import etomica.virial.GUI.views.ParameterWindowView;





public class ParameterController {
	
	private ParameterWindowView PView;
	private SpeciesModel SModel;
	private SingleSpeciesModel SSModel;
	private LJDefaultParametersModel LJDModel;
	
	public ParameterController(ParameterWindowView pView,SuperModel SM){
		PView = pView;
		SModel = SM.getSm();
		SSModel = SM.getSsm();
		LJDModel = SM.getLJD();
	/*	
		PView.getSpeciesParameters().addChoosePureRadio(new ChoosePureRButtonListener());
		PView.getSpeciesParameters().addChooseBinaryRadio(new ChooseBinaryRButtonListener());
		PView.getOneSpecies().addChooseComponentListener(new ChooseComponentListener());
		PView.getOneSpecies().addChooseLJListener(new ChooseLJListener());
		PView.getOneSpecies().addChooseCO2Listener(new ChooseCO2Listener());
		PView.getOneSpecies().addChooseResetListener(new ChooseResetListener());
		PView.getOneSpecies().addChooseDefaultListener(new ChooseDefaultListener());*/
		//PView.getDialogBox().addChooseOkayButtonListener(new ChooseOkayButtonListener());
	//	PView.getTwoSpecies().addChooseComponentAListener(new ChooseComponentAListener());
		
		
	}
/*	
	class ChoosePureRButtonListener implements ActionListener{
		public void actionPerformed(ActionEvent E){
			System.out.println("You ve chosen To run a simulation with Pure Species");
			CardLayout cl = (CardLayout)PView.getSpeciesParameters().getCards().getLayout();
			cl.show(PView.getSpeciesParameters().getCards(), "PURE");
		}
		
	}
	class ChooseBinaryRButtonListener implements ActionListener{
		public void actionPerformed(ActionEvent E){
			System.out.println("You ve chosen To run a simulation with Binary Mixture");
			CardLayout cl = (CardLayout)PView.getSpeciesParameters().getCards().getLayout();
			cl.show(PView.getSpeciesParameters().getCards(), "BINARY");

		}
		
	}
	
	class ChooseComponentListener implements ActionListener{
		public void actionPerformed(ActionEvent E){
			 JComboBox cb = (JComboBox)E.getSource();
			 String SelectedComponent = (String)cb.getSelectedItem();
			 System.out.println(SelectedComponent);
			 System.out.println(PView.getOneSpecies());
			 CardLayout cl = (CardLayout)PView.getOneSpecies().getPotentialCards().getLayout();
			 PView.getOneSpecies().getPotentialCards().setVisible(true);
			 SingleSpeciesParam.getPotentialLabel().setVisible(true);
			 if(SelectedComponent == "LJ"){
	         cl.show(PView.getOneSpecies().getPotentialCards(), "LJ");
	         
			 }
			 
			 if(SelectedComponent == "CO2"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "CO2");}
			
			 if(SelectedComponent == "Methanol"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "Alcohol");}
			
			 if(SelectedComponent == "Ethanol"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "Alcohol");}
				
			 
			 if(SelectedComponent == "Methane"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
			
			 if(SelectedComponent == "Ethane"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
				
			 if(SelectedComponent == "Propane"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
				
		 }
	}
	
	class ChooseLJListener implements ActionListener{
		public void actionPerformed(ActionEvent E){
			JComboBox cb = (JComboBox)E.getSource();
			String SelectedComponent = (String)cb.getSelectedItem();
			System.out.println(SelectedComponent);
				PView.getOneSpecies().setSigmaHSRefLJText(SSModel.getSigmaHSRefLJ());
				PView.getOneSpecies().getSigmaHSRefLJField().setVisible(true);
				SingleSpeciesParam.getSigmaHSRefLJLabel().setVisible(true);
				PView.getOneSpecies().getButtonText().setVisible(true);
				PView.getOneSpecies().getDefaultValues().setVisible(true);
				PView.getOneSpecies().getResetLabel().setVisible(true);
				PView.getOneSpecies().getReset().setVisible(true);
		}
	}
	
	class ChooseCO2Listener implements ActionListener{
		public void actionPerformed(ActionEvent E){
			JComboBox cb = (JComboBox)E.getSource();
			String SelectedComponent = (String)cb.getSelectedItem();
			System.out.println(SelectedComponent);
				if(SelectedComponent == "2-Centered-with-Quad"){
					PView.getOneSpecies().setSigmaHSRefLJText(SSModel.getSigmaHSRefCO2LJ());
					PView.getOneSpecies().getSigmaHSRefLJField().setVisible(true);
					SingleSpeciesParam.getSigmaHSRefLJLabel().setVisible(true);
					PView.getOneSpecies().getButtonText().setVisible(true);
					PView.getOneSpecies().getDefaultValues().setVisible(true);
					PView.getOneSpecies().getResetLabel().setVisible(true);
					PView.getOneSpecies().getReset().setVisible(true);
				}
				
				if(SelectedComponent == "EPM2-with-LJ"){
					PView.getOneSpecies().setSigmaHSRefLJText(SSModel.getSigmaHSRefCO2EPM2());
					PView.getOneSpecies().getSigmaHSRefLJField().setVisible(true);
					SingleSpeciesParam.getSigmaHSRefLJLabel().setVisible(true);
					PView.getOneSpecies().getButtonText().setVisible(true);
					PView.getOneSpecies().getDefaultValues().setVisible(true);
					PView.getOneSpecies().getResetLabel().setVisible(true);
					PView.getOneSpecies().getReset().setVisible(true);
				}
				
				if(SelectedComponent == "TRAPPE-with-LJ"){
					//default title and icon
					PView.getDialogFrame().setVisible(true);
					PView.getDialogFrame().setBounds(300,300, 300,100);
				}
		}
	}
		
		class ChooseOkayButtonListener implements ActionListener{
			public void actionPerformed(ActionEvent E){
				PView.getDialogFrame().setVisible(false);
				PView.reset();
				//PView.getDialogFrame().dispose();
			}
			
		}
		
		
		class ChooseResetListener implements ActionListener{
			
			public void actionPerformed(ActionEvent E){
				PView.reset();
				//PView.getDialogFrame().dispose();
			}
			
		} 
		
		class ChooseDefaultListener implements ActionListener{
			
			public void actionPerformed(ActionEvent E){
				int SelectedComponent = PView.getOneSpecies().getComponentsList().getSelectedIndex();
				
				if(SelectedComponent == 0)
				{
					int SelectedPotential = PView.getOneSpecies().getLJpotentialList().getSelectedIndex();
					if(SelectedPotential == 0){
						PView.getDefaultViewFrame().setVisible(true);
						PView.getDefaultViewFrame().setBounds(300,300, 300,400);
						PView.getLJDefault().getSigmaField1().setText(Double.toString(LJDModel.getSigmaValue()));
						PView.getLJDefault().getEpsilonField1().setText(Double.toString(LJDModel.getEpsilonValue()));
						CardLayout cl = (CardLayout)(PView.getLJDefault().getLJDefaultViewCards().getLayout());
						cl.show(PView.getLJDefault().getLJDefaultViewCards(), "LJ1");
					}
					if(SelectedPotential == 1){
						PView.getDefaultViewFrame().setVisible(true);
						PView.getDefaultViewFrame().setBounds(300,300, 300,400);
						PView.getLJDefault().getSigmaField2().setText(Double.toString(LJDModel.getSigmaValue()));
						PView.getLJDefault().getEpsilonField2().setText(Double.toString(LJDModel.getEpsilonValue()));
						PView.getLJDefault().getMomentField().setText(Double.toString(LJDModel.getMomentValue()));
						PView.getLJDefault().getBondLField().setText(Double.toString(LJDModel.getBondLenthValue()));
						CardLayout cl = (CardLayout)(PView.getLJDefault().getLJDefaultViewCards().getLayout());
						cl.show(PView.getLJDefault().getLJDefaultViewCards(), "LJ2");
					}
					
				}
				if(SelectedComponent == 1)
				{
					int SelectedPotential = PView.getOneSpecies().getCO2PotentialList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				if(SelectedComponent == 2)
				{
					int SelectedPotential = PView.getOneSpecies().getAlcoholPotentialsList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				if(SelectedComponent == 3)
				{
					int SelectedPotential = PView.getOneSpecies().getAlcoholPotentialsList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				
				if(SelectedComponent == 4)
				{
					int SelectedPotential = PView.getOneSpecies().getAlkanePotentialList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				if(SelectedComponent == 5)
				{
					int SelectedPotential = PView.getOneSpecies().getAlkanePotentialList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				if(SelectedComponent == 6)
				{
					int SelectedPotential = PView.getOneSpecies().getAlkanePotentialList().getSelectedIndex();
					System.out.println(SelectedPotential);
				}
				
			}
			
		} 
		
		class ChooseComponentAListener implements ActionListener{
			public void actionPerformed(ActionEvent E){
				 JComboBox cb = (JComboBox)E.getSource();
				 String SelectedComponent = (String)cb.getSelectedItem();
				 System.out.println(SelectedComponent);
				 System.out.println(PView.getTwoSpecies());
				 CardLayout cl = (CardLayout)PView.getTwoSpecies().getComponentBCards().getLayout();
				 PView.getTwoSpecies().getComponentBCards().setVisible(true);
				 if(SelectedComponent == "LJ"){
		         cl.show(PView.getOneSpecies().getPotentialCards(), "CompBLJ");
		         
				 }
				 
				 if(SelectedComponent == "CO2"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "CO2");}
				
				 if(SelectedComponent == "Methanol"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "Alcohol");}
				
				 if(SelectedComponent == "Ethanol"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "Alcohol");}
					
				 
				 if(SelectedComponent == "Methane"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
				
				 if(SelectedComponent == "Ethane"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
					
				 if(SelectedComponent == "Propane"){
			         cl.show(PView.getOneSpecies().getPotentialCards(), "Alkane");}
					
			 }
		}
		
		*/
		
	
	
}
