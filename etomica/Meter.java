package simulate;

import java.awt.*;
import java.util.Vector;
import java.util.*;

public abstract class Meter extends Component
{
    int updateInterval = 1;
    int iieCount = 0;
    boolean firstCall = true;
    double sum, sumSquare;
    int count = 0;
    private Meter nextMeter, previousMeter;
    public DataDisplay parentDisplay;

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
    
	public void update(Phase phase) {
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
	    return (count>0) ? sum/(double)count : Double.MAX_VALUE;
	}
	
	public double error() {
	    double avg = average();
	    return (count>1) ? Math.sqrt((sumSquare/(double)count - avg*avg)/(double)(count-1)) : Double.MAX_VALUE;
	}
	
	public void reset() {
	    count = 0;
	    sum = 0;
	    sumSquare = 0;
	}
	
	public abstract double currentValue(Phase phase);

    public final void setNextMeter(Meter meter) {
      this.nextMeter = meter;
      if(meter != null) {meter.previousMeter = this;}
    }
    
    public final Meter getNextMeter() {return nextMeter;}
    public final Meter getPreviousMeter() {return previousMeter;}

}	 
