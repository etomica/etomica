package simulate;

import java.awt.Component;

public abstract class Meter extends Component implements IntegrationIntervalListener
{
    int updateInterval;
    int iieCount;
    boolean firstCall = true;
    double sum, sumSquare;
    int count = 0;
    private Meter nextMeter, previousMeter;
    public Phase phase;
    String label;

	public Meter() {
	    setUpdateInterval(1);
	    label = "Property";
	}
	
	public String getLabel() {
	    return label;
	}
	public void setLabel(String s) {
	    label = s;
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
	        double value = currentValue();
	        if(!Double.isNaN(value)) {
	            sum += value;
	            sumSquare += value*value;
	            count++;
	        }
	    }
	}
		
	public double average() {
	    return (count>0) ? sum/(double)count : 0.0;
	}
	
	public double error() {    //temporary---needs to be rewritten to do block averaging
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
    
    public void initialize() {;}
    
    public final Meter getNextMeter() {return nextMeter;}
    public final Meter getPreviousMeter() {return previousMeter;}
    
    public void integrationIntervalAction(IntegrationIntervalEvent evt) {updateStatistics(phase);}

}	 
