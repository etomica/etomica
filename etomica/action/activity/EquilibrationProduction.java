/*
 * History
 * Created on Oct 29, 2004 by kofke
 */
package etomica.action.activity;

import etomica.Action;
import etomica.ActivityIntegrate;
import etomica.Integrator;
import etomica.action.ResetAccumulators;
import etomica.utility.java2.LinkedList;


/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class EquilibrationProduction extends ActivityGroupSeries {

	public EquilibrationProduction(Integrator equilibrationIntegrator,
			Integrator productionIntegrator, LinkedList accumulatorManagerList) {
		equilibrationActivity = new ActivityIntegrate(equilibrationIntegrator);
		productionActivity = new ActivityIntegrate(productionIntegrator);
		equilibrationActivity.getIntegrator().setEquilibrating(true);
		addAction(equilibrationActivity);
		
		productionPreparationActivity = new ActivityGroupSeries();
		productionPreparationActivity.addAction(new ResetAccumulators(accumulatorManagerList));
		productionPreparationActivity.addAction(setIntegratorForProduction);
		addAction(productionPreparationActivity);
		
		addAction(productionActivity);
	}
	
	public EquilibrationProduction(Integrator commonIntegrator, LinkedList accumulatorManagerList) {
		this(commonIntegrator,commonIntegrator,accumulatorManagerList);
	}

	public ActivityIntegrate getEquilibrationActivity() {
		return equilibrationActivity;
	}
	public ActivityGroupSeries getProductionPreparationActivity() {
		return productionPreparationActivity;
	}
	public ActivityIntegrate getProductionActivity() {
		return productionActivity;
	}
	
	private final ActivityIntegrate equilibrationActivity, productionActivity;
	private final ActivityGroupSeries productionPreparationActivity;
	
	private Action setIntegratorForProduction = new Action() {
		private String label = "Set integrator for production";
		public void actionPerformed() {
			productionActivity.getIntegrator().setEquilibrating(false);
		}
		public String getLabel() {return label;}
		public void setLabel(String label) {this.label = label;}
	};
}
