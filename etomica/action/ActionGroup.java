package etomica.action;

public interface ActionGroup {

    public boolean removeAction(Action oldAction);
    
    public void addAction(Action newAction);
    
	public Action[] getAllActions();

}