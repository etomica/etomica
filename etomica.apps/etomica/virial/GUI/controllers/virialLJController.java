package etomica.virial.GUI.controllers;

import etomica.virial.GUI.models.virialLJModel;
import etomica.virial.GUI.views.virialLJModelView;
import etomica.virial.simulations.VirialLJ;
import etomica.virial.simulations.VirialLJ.VirialLJParam;
import java.awt.event.*;
import java.io.*;
public class virialLJController {

	 //Instances of Model and View
	 private virialLJModel lj_model;
	 private virialLJModelView  lj_view;

	 //Constructor
	 public virialLJController(virialLJModel model, virialLJModelView view) {
		 	lj_model = model;
		 	lj_view  = view;
	        
	        //listeners to the view.
	        view.addSendParametersListener(new SendParametersListener());
	        view.addResetParametersListener(new ResetParametersListener());
	    }
	 
	//Listener class to run the VirialLJ simulation when the 'Run' button is pressed
	 class SendParametersListener implements ActionListener 
	 {
	        public void actionPerformed(ActionEvent e) 
	        {
	        	int nPoints;
	            double temperature;
	            long steps;
	            double sigmaHSRef;
	            
	            try {
	            	nPoints = lj_view.getUserInputNPoints();
	            	temperature = lj_view.getUserInputTemperature();
	            	steps = lj_view.getUserInputSteps();
	            	sigmaHSRef = lj_view.getUserInputSigmaHSRef();
	            	Writer output = null;
	            	String nPointsText = "nPoints " + nPoints;
	            	String temperatureText = "temperature " + temperature;
	            	String stepsText = "numSteps " + steps;
	            	String sigmaHSRefText = "sigmaHSRef " + sigmaHSRef;
	            	File file = new File("input.txt");
	            	//VirialLJ tempVirialLJ = new VirialLJ();
	            	//VirialLJ.runVirial(new VirialLJParam(nPoints,temperature,steps,sigmaHSRef));
	            	
	            	try{
	                output = new BufferedWriter(new FileWriter(file));
	                output.write(nPointsText);
	                output.write("\n");
	                output.write(temperatureText);
	                output.write("\n");
	                output.write(stepsText);
	                output.write("\n");
	                output.write(nPointsText);
	                output.write("\n");
	                output.write(sigmaHSRefText);
	                output.close();
	                System.out.println("Your file has been written"); }
	            	catch (Exception E){
	            		
	            	}
	            	
	                //lj_model.sendParameters(userInput);
	                
	            } catch (NumberFormatException nfex) {
	            	lj_view.showError("Bad input:");
	            }
	        }
	    }//end inner class SendParametersListener

	 class ResetParametersListener implements ActionListener {
	        public void actionPerformed(ActionEvent e) {
	        	lj_model.reset();
	        	lj_view.reset();
	        }
	    }// end inner class ResetParametersListener

}
