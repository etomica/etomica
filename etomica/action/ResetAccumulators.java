/*
 * History
 * Created on Nov 4, 2004 by kofke
 */
package etomica.action;

import etomica.AccumulatorManager;
import etomica.Action;
import etomica.utility.java2.Iterator;
import etomica.utility.java2.LinkedList;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
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
			((AccumulatorManager)iterator.next()).resetAccumulators();
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
