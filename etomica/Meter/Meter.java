package simulate;

import java.awt.*;
import java.util.Vector;
import java.util.*;

public abstract class Meter implements IntegrationIntervalListener
{
    int updateInterval = 1;
    int iieCount = 0;
    boolean firstCall = true;

	public Meter()
	{
	}

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(updateInterval > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
	public void update(IntegrationIntervalEvent iie) {
	    if(firstCall) {
	        firstCall = false;
	        initialEnergy = iie.phase.getTotalEnergy();
	    }
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        takeMeasurement(iie.phase);
	        double totalEnergy = phase.getTotalEnergy();
	        double error = 100*(totalEnergy-initialEnergy)/initialEnergy;
	        System.out.println(error);
	    }
	}
	
	public abstract void takeMeasurement(Phase phase);
}	 
