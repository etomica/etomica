/*
 * History
 * Created on Oct 29, 2004 by kofke
 */
package etomica.action.activity;

import etomica.ActivityIntegrate;
import etomica.Integrator;


/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class EquilibrationProduction extends ActivityGroupSeries {

	public EquilibrationProduction(Integrator equilibrationIntegrator,
			Integrator productionIntegrator) {
		equilibrationActivity = new ActivityIntegrate(equilibrationIntegrator);
		productionActivity = new ActivityIntegrate(productionIntegrator);
		addAction(equilibrationActivity);
		addAction(productionActivity);
	}
	
	public EquilibrationProduction(Integrator commonIntegrator) {
		this(commonIntegrator,commonIntegrator);
	}

	public ActivityIntegrate getEquilibrationActivity() {
		return equilibrationActivity;
	}
	public ActivityIntegrate getProductionActivity() {
		return productionActivity;
	}
	
	private final ActivityIntegrate equilibrationActivity, productionActivity;
}
