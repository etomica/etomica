package etomica.listener;

import etomica.action.IAction;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;

public class IntegratorListenerAction implements IIntegratorListener {

    private IAction action;
    private int interval;
    private int intervalCount;
    
    public IntegratorListenerAction(IAction a) {
        this(a, 1);
    }
    
    public IntegratorListenerAction(IAction a, int interval) {
        action = a;
        this.interval = interval;
        intervalCount = 0;
    }
    
    public void integratorInitialized(IIntegratorEvent e) {
        if(intervalCount >= interval) {
            
        }
    }
    
    public void integratorStepStarted(IIntegratorEvent e) {
        if(++intervalCount >= interval) {
            
        }
        
    }
    
    public void integratorStepFinished(IIntegratorEvent e) {
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
