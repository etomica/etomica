package etomica.action;
import etomica.MeterAbstract;
import etomica.Simulation;
import etomica.utility.java2.Iterator;
import etomica.utility.java2.LinkedList;

 /**
  * Superclass of classes that apply some elementary action (transformation) to 
  * one or more Meters.
  *
  * @author David Kofke
  * 
  */
 
 /* History
 *  12/31/02  (DAK) changed package to etomica.action
 */
public abstract class MeterAction implements etomica.Action {

    protected MeterAbstract[] meters;
    protected LinkedList meterList;
    protected String label = "Meter Action";
    
    /**
     * Constructor giving array defining the meters subject to this action.
     */
    public MeterAction(MeterAbstract[] m) {
        super();
        meters = m;
    }
    /**
     * Meterlist of given simulation defines meters subject to this action.
     */
    public MeterAction(Simulation sim) {
        super();
        meterList = sim.getMeterList();
    }
        
    /**
     * Returns the array of all meters subject to action.
     */
    public MeterAbstract[] getMeters() {return meters;}
    
    /**
     * Returns a descriptive label for this action.
     */
	public String getLabel() {
		return label;
	}
	/**
	 * Sets a descriptive label for this action.  This might be referenced,
	 * for example, by a button invoking this action in a graphical interface.
	 */
	public void setLabel(String label) {
		this.label = label;
	}
    /**
     * Sets the ith meter in the array of meters.  Takes
     * no action if i is greater than the current number of meters (minus 1).
     */
    public MeterAbstract getMeters(int i) {
        if(meters == null) return null;
        if(i < 0 || i >= meters.length) throw new ArrayIndexOutOfBoundsException();
        else return meters[i];
    }
    /**
     * Sets the array of meters being tabulated.  Existing sources
     * are discarded.
     */
    public void setMeters(MeterAbstract[] s) {
        if(s == null || s.length == 0 || s[0] == null) {
            meters = null;
            return;
        }
        meters = s;
    }
    
    public void setMeters(int i, MeterAbstract s) {
        if(meters == null && i != 0) throw new NullPointerException();
        else if(meters == null && i == 0) {setMeters(s); return;}
        if(i < 0 || i >= meters.length) throw new ArrayIndexOutOfBoundsException();
        else meters[i] = s;
    }
    
    /**
     * Sets the given meter as the only meter subject to action.
     * Existing sources are discarded.
     */
    public void setMeters(MeterAbstract s) {
        setMeters(new MeterAbstract[] {s});
    }
    
    
    /**
     * Adds the given source to the sources being plotted.
     * Existing sources are retained.
     */
    public void addMeters(MeterAbstract s) {
        if(s == null) return;
        int nMeter = (meters == null) ? 0 : meters.length;
        MeterAbstract[] newMeters = new MeterAbstract[nMeter+1];
        for(int i=0; i<nMeter; i++) newMeters[i] = meters[i];
        newMeters[nMeter] = s;
        setMeters(newMeters);
    }
    
    /**
     * Adds the given sources to the sources being plotted.
     * Existing sources are retained.
     */
    public void addMeters(MeterAbstract[] s) {
        if(s == null) return;
        int nMeter = (meters == null) ? 0 : meters.length;
        MeterAbstract[] newMeters = new MeterAbstract[nMeter+s.length];
        for(int i=0; i<nMeter; i++) newMeters[i] = meters[i];
        for(int i=nMeter; i<newMeters.length; i++) newMeters[i] = s[i-nMeter];
        setMeters(newMeters);
    }
        
    /**
     *  Performs action on all meters.
     */
    public void actionPerformed() {
        if(meters != null) {
            for(int i=0; i<meters.length; i++) {
                actionPerformed(meters[i]);
            }
        }
        if(meterList != null) {
            for(Iterator iter=meterList.iterator(); iter.hasNext(); ) {
                actionPerformed((MeterAbstract)iter.next());
            }
        }
    }
    
    /**
     * Defines action to be performed on each meter.
     */
    public abstract void actionPerformed(MeterAbstract meter);

        

}