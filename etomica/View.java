package simulate;

import java.awt.*;
import simulate.*;
import java.beans.*;

public abstract class View extends java.awt.Component
{
    int updateInterval;
    int iieCount;
    public Display parentDisplay;
    
	public View() {
	    setUpdateInterval(1);
	}
	
    public abstract void doUpdate();
    
    public final void updateView() {
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
	    }
	}   

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
}