package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;

public class DataDisplay extends Display
{
    
	public DataDisplay()
	{

		//{{INIT_CONTROLS
		setLayout(null);
		Insets ins = getInsets();
		setSize(ins.left + ins.right + 430,ins.top + ins.bottom + 270);
		//}}
	}

	//{{DECLARE_CONTROLS
	Meter firstMeter, lastMeter;
	private int nMeters = 0;
	//}}
	
	public void setPhase(Phase p) {phase = p;}
	
	public void add(Meter m) {
	    if(firstMeter == null) {firstMeter = m;}
	    if(lastMeter != null) {lastMeter.setNextMeter(m);}
	    lastMeter = m;
	    nMeters++;
	    m.parentDisplay = this;
	}
	
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= nMeters) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.getNextMeter();}  //get ith meter in list
        return m;
    }

    public void updateData(IntegrationIntervalEvent evt)
    {
        // This method is derived from interface simulate.IntegrationIntervalListener
        for(Meter m=firstMeter; m!=null; m=m.getNextMeter()) {
            m.updateStatistics(evt.phase);
        }
        view.updateView();
        repaint();
    }
}