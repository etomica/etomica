package simulate;

import java.awt.*;
import java.util.Vector;
import java.util.*;

public abstract class Meter extends Component implements IntegrationIntervalListener
{
    int updateInterval;
    int iieCount;
    boolean firstCall = true;
    double sum, sumSquare;
    int count = 0;
    private Meter nextMeter, previousMeter;
    public DataDisplay parentDisplay;
    public Phase phase;

	public Meter() {
	    setUpdateInterval(1);
	}

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
	public void updateStatistics(Phase phase) {
	    if(firstCall) {
	        firstCall = false;
	    }
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        double value = currentValue(phase);
	        sum += value;
	        sumSquare += value*value;
	        count++;
	    }
	}
		
	public double average() {
	    return (count>0) ? sum/(double)count : 0.0;
	}
	
	public double error() {
	    double avg = average();
	    return (count>1) ? Math.sqrt((sumSquare/(double)count - avg*avg)/(double)(count-1)) : 0.0;
	}
	
	public void reset() {
	    count = 0;
	    sum = 0;
	    sumSquare = 0;
	}
	
	public abstract double currentValue();

    public final void setNextMeter(Meter meter) {
      this.nextMeter = meter;
      if(meter != null) {meter.previousMeter = this;}
    }
    
    public final Meter getNextMeter() {return nextMeter;}
    public final Meter getPreviousMeter() {return previousMeter;}
    
    public void updateData(IntegrationIntervalEvent evt) {updateStatistics(phase);}

}	 
