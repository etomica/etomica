/*
 * History
 * Created on Nov 4, 2004 by kofke
 */
package etomica.action;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.Action;
import etomica.data.DataAccumulator;

/**
 * Action that performs a call to the reset() method of a set
 * of accumulators, as specified via a list of AccumulatorManager
 * instances.
 */
public class ResetAccumulators implements Action {

	/**
	 * 
	 */
	public ResetAccumulators(LinkedList accumulatorManagerList) {
		this.accumulatorManagerList = accumulatorManagerList;
	}

	public void actionPerformed() {
		Iterator iterator = accumulatorManagerList.iterator();
		while (iterator.hasNext()) {
			((DataAccumulator)iterator.next()).reset();
		}
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	private final LinkedList accumulatorManagerList;
	private String label = "Reset Accumulators";
}
