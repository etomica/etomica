/*
 * History
 * Created on Oct 29, 2004 by kofke
 */
package etomica.action.activity;

import java.util.LinkedList;

import etomica.action.Action;
import etomica.action.ResetAccumulators;
import etomica.action.SimulationDataAction;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorPhase;
import etomica.simulation.Simulation;

/**
 * Pre-configured activity group that runs the integrator through a set of
 * equilibration cycles, then sets integrator to a non-equilibrating mode, and
 * runs integrator through a set of production cycles. Between equilibration and
 * production, activity will reset a given list of accumulators, which typically
 * would be obtained by passing the result of a call to
 * simulation.getAccumulatorManagerList.
 */
public class EquilibrationProduction extends ActivityGroupSeries {

	/**
	 * Constructs activity such that a two different integrators instances are
	 * used for the equilibration and production integrations.
	 * 
	 * @param equilibrationIntegrator
	 *            the integrator used for equilibration
	 * @param productionIntegrator
	 *            the integrator used for production
	 * @param accumulatorManagerList
	 *            list of accumulatorManagers that will be reset between
	 *            equilibration and production periods
	 */
	public EquilibrationProduction(Simulation sim, IntegratorPhase equilibrationIntegrator,
			IntegratorPhase productionIntegrator, LinkedList dataManagerList) {
		equilibrationActivity = new ActivityIntegrate(sim,equilibrationIntegrator);
		productionActivity = new ActivityIntegrate(sim,productionIntegrator);
		equilibrationIntegrator.setEquilibrating(true);
		addAction(equilibrationActivity);

		productionPreparationActivity = new ActivityGroupSeries();
		productionPreparationActivity.addAction(new SimulationDataAction(sim,
				new ResetAccumulators()));
		productionPreparationActivity.addAction(new MyAction(productionIntegrator));
		addAction(productionPreparationActivity);

		addAction(productionActivity);
	}

	/**
	 * Constructs activity such that a common integrator instance is used for
	 * both the equilibration and production integrations.
	 * 
	 * @param commonIntegrator
	 *            the integrator used for equilibration and production
	 * @param accumulatorManagerList
	 *            list of accumulatorManagers that will be reset between
	 *            equilibration and production periods
	 */
	public EquilibrationProduction(Simulation sim, IntegratorPhase commonIntegrator,
			LinkedList accumulatorManagerList) {
		this(sim, commonIntegrator, commonIntegrator, accumulatorManagerList);
	}

	/**
	 * @return the activity instance that performs the equilibration steps
	 */
	public ActivityIntegrate getEquilibrationActivity() {
		return equilibrationActivity;
	}

	/**
	 * @return the group of activities that are performed between the
	 *         equilibration and production periods
	 */

	public ActivityGroupSeries getProductionPreparationActivity() {
		return productionPreparationActivity;
	}

	/**
	 * @return the activity instance that performs the production steps
	 */
	public ActivityIntegrate getProductionActivity() {
		return productionActivity;
	}

	/**
	 * Sets the number of integrator steps to be performed during the
	 * equilibration period.
	 */
	public void setEquilibrationSteps(int nSteps) {
		equilibrationActivity.setMaxSteps(nSteps);
	}

	/**
	 * @return the number of steps integrator performed during the production
	 *         period.
	 */
	public long getProductionSteps() {
		return productionActivity.getMaxSteps();
	}

	/**
	 * Sets the number of integrator steps to be performed during the production
	 * period.
	 */
	public void setProductionSteps(int nSteps) {
		productionActivity.setMaxSteps(nSteps);
	}

	/**
	 * @return the number of integrator steps performed during the equilibration
	 *         period.
	 */
	public long getEquilibrationSteps() {
		return equilibrationActivity.getMaxSteps();
	}

	private final ActivityIntegrate equilibrationActivity, productionActivity;

	private final ActivityGroupSeries productionPreparationActivity;

	/**
	 * Action that is part of group between equilibration and production, and
	 * which results in the integrators being set for production.  This cannot
	 * be done at construction because equilibriation and production might use the
	 * same integrator instance.
	 */
	private static class MyAction implements Action {
		private String label = "Set integrator for production";
        private Integrator integrator;

        public MyAction(Integrator integrator) {
            this.integrator = integrator;
        }
        
		public void actionPerformed() {
			integrator.setEquilibrating(false);
		}

		public String getLabel() {
			return label;
		}

		public void setLabel(String label) {
			this.label = label;
		}
	};
}