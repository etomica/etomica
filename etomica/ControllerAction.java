package etomica;

 /**
  * Elementary action performed on a controller.
  */
public abstract class ControllerAction implements etomica.Action, ControllerListener {

	private String label = "ControllerAction";
	protected Controller controller;
    public ControllerAction() {this(null);}
    public ControllerAction(Controller c) {
        super();
        controller = c;
    }
        
	public String getLabel() {
		return label;
	}
	public void setLabel(String label) {
		this.label = label;
	}
    public void setController(Controller c) {controller = c;}
    public Controller getController() {return controller;}
    
    public void actionPerformed(SimulationEvent evt) {
        actionPerformed((ControllerEvent)evt);
    }
    public void actionPerformed(ControllerEvent ce) {
        actionPerformed(ce.controller());
    }
    
    public void actionPerformed() {actionPerformed(controller);}
    
    public abstract void actionPerformed(Controller c);

}