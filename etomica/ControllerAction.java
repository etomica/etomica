package etomica;

 /**
  * Elementary action performed on a controller.
  */
public abstract class ControllerAction extends etomica.Action implements ControllerListener {

    public static String getVersion() {return "ControllerAction:01.11.20/"+Action.VERSION;}

    protected Controller controller;
    public ControllerAction() {this(null);}
    public ControllerAction(Controller c) {
        super();
        controller = c;
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