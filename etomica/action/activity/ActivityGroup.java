package etomica.action.activity;

import etomica.action.Action;
import etomica.action.ActionGroup;

public interface ActivityGroup extends ActionGroup {

	public Action[] getCompletedActions();
	public Action[] getCurrentActions();
	public Action[] getPendingActions();
}
