package etomica;

public class EtomicaInfo {
    
    private String description;
    private boolean enabled = true;
    
    public EtomicaInfo() {
        this("No description available");
    }
    public EtomicaInfo(String desc) {
        description = desc;
    }
    
    public String getDescription() {return description;}
    public void setDescription(String str) {description = str;}
    
    /**
     * Accessor for flag indicating if component is enabled for use in Etomica.  May be set
     * to false if other required classes are absent (for example, 3D graphics libraries
     * needed for display).  Default is <code>true</code>.
     */
    public boolean isEnabled() {return enabled;}
    /**
     * Mutator for flag indicating if component is enabled for use in Etomica.  May be set
     * to false if other required classes are absent (for example, 3D graphics libraries
     * needed for display).  Default is <code>true</code>.
     */
    public void setEnabled(boolean b) {enabled = b;}
}