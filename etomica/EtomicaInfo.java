package etomica;

import java.util.HashMap;
import java.util.LinkedList;

import etomica.compatibility.FeatureSet;


public class EtomicaInfo {
    
    private String description;
    private boolean enabled = true;
    
    public EtomicaInfo() {
        this("No description available");
    }
    public EtomicaInfo(String desc) {
        description = desc;
    }
   
	public FeatureSet getFeatures()
	{
		FeatureSet fts = new FeatureSet();
		return fts.add( "CLASS", this.getClass().getName() )
		.add( "DESCRIPTION", description )
		.add( "SHORT_DESCRIPTION", getShortDescription() );
	}
    protected static final int SHORT_DESCRIPTION_LENGTH = 50;
    
    /** The need for a short text to put on pull down menus made it necessary to add this function. By default
     * it takes the description field (possibly long) and will return a substring of it, indicating continuation with
     * with ellipsis (...)
     * @return Short description of this field
     */ 
    public String getShortDescription() 
    { 
    	if ( description.length()<SHORT_DESCRIPTION_LENGTH )
    		return description;
    	return description.substring( SHORT_DESCRIPTION_LENGTH );
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