package etomica.listener;

import etomica.api.IAction;
import etomica.api.IEvent;
import etomica.api.IIntegratorListener;

public class IntegratorListenerAction implements IIntegratorListener {

    private IAction action;
    private int interval;
    private int intervalCount;
    
    public IntegratorListenerAction(IAction a) {
        action = a;
        interval = 1;
        intervalCount = 0;
    }
    
    public void integratorInitialized(IEvent e) {
        if(intervalCount >= interval) {
            
        }
    }
    
    public void integratorStepStarted(IEvent e) {
        if(++intervalCount >= interval) {
            
        }
        
    }
    
    public void integratorStepFinished(IEvent e) {
        if(intervalCount >= interval) {
            intervalCount = 0;
            action.actionPerformed();
        }
    }
    
    public void setInterval(int i) {
        interval = i;
    }
    
    public int getInterval() {
        return interval;
    }
    
    public IAction getAction() {
        return action;
    }
    
}
