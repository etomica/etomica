package etomica.action;
import etomica.Meter;

 /**
  * Interface for object that performs action on a specified Meter.
  *
  * @author David Kofke
  * 
  */
 
public interface MeterAction extends etomica.Action {

    /**
     * Returns the meter subject to action.
     */
    public Meter getMeter();
    
    /**
     * Sets the given meter as the one subject to action.
     */
    public void setMeter(Meter s);

}