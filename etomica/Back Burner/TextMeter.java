package simulate;

import java.awt.*;
import java.util.Vector;
import java.util.*;

public class TextMeter extends java.awt.Panel implements IntegrationIntervalListener
{
    int updateInterval = 1;
    int iieCount = 0;
    boolean firstCall = true;
    double initialEnergy;

	public TextMeter()
	{
		setLayout(null);
		setSize(200,200);
	}

    public int getUpdateInterval() {return updateInterval;}
    public void setUpdateInterval(int i) {updateInterval = i;}
    
	public void updateAverage(IntegrationIntervalEvent iie) {
	    if(firstCall) {
	        firstCall = false;
	        initialEnergy = iie.phase.getTotalEnergy();
	    }
	    iieCount++;
	    if(iieCount == updateInterval) {
	        iieCount = 0;
	        Phase phase = iie.phase;
	        double totalEnergy = phase.getTotalEnergy();
	        double error = 100*(totalEnergy-initialEnergy)/initialEnergy;
	        System.out.println(error);
	    }
	}
}	 
