package etomica;

/**
 * Marker interface to indicate that a class is suitable to be displayed as an
 * option for addition to a simulation via the Etomica GUI.  Class must implement
 * this interface to be given as a radio-button choice in one of the simulation editor panes.
 */
public interface EtomicaElement {
    /**
     * Method to set the name of this Meter
     * 
     * @param name The name string to be associated with this Meter
     */
    public void setName(String name);

    /**
     * Accessor method of the name of this Meter
     * 
     * @return The given name of this Meter
     */
    public String getName();

}
