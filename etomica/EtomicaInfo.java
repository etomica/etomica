package etomica;

import java.lang.reflect.InvocationTargetException;

import etomica.compatibility.FeatureSet;


public class EtomicaInfo implements java.io.Serializable {
    
    protected static final int MAX_SHORT_DESCRIPTION_LENGTH = 50;
   
    private String description;
    private String short_description;
    private boolean enabled = true;
    private FeatureSet features;
    
    public EtomicaInfo() {
        this("No description available");
        setDefaultFeatureSet();
    }
    
    protected void setDefaultFeatureSet()
    {
    	features = new FeatureSet();
    	features.add( "CLASS", this.getClass().getName() )
		        .add( "DESCRIPTION", description )
		        .add( "SHORT_DESCRIPTION", short_description );
    }

	public EtomicaInfo(String desc) {
        description = desc;
        if ( desc.length()>MAX_SHORT_DESCRIPTION_LENGTH ) 
        	short_description = desc.substring(0, MAX_SHORT_DESCRIPTION_LENGTH-3)+"..."; 
        else
        	short_description = desc;
        setDefaultFeatureSet();
    }
   
    public EtomicaInfo(String desc, String short_desc ) {
        description = desc;
        short_description = short_desc;
        setDefaultFeatureSet();
    }

    public FeatureSet getFeatures()
	{
    	return features;
	}
	
	public static EtomicaInfo getInfo( Class myclass )
	{
		etomica.EtomicaInfo info = null;
        java.lang.reflect.Method method = null;
        
        try {
            method = myclass.getMethod("getEtomicaInfo",null);
        }
        catch ( NoSuchMethodException se ) {
        	// class does not have a getEtomicaInfo method
        	return new EtomicaInfo( myclass.getName() );
        }

        try {
            info = (etomica.EtomicaInfo)method.invoke(myclass, null);
        }
        catch (IllegalAccessException e) {
            System.out.println("Exception retrieving info for class " + myclass.getName()+ ": " + e.getLocalizedMessage() );
        }
        catch (InvocationTargetException e) {
            System.out.println("Exception retrieving info for class " + myclass.getName()+ ": " + e.getLocalizedMessage() );
        }
        
        return info;
	}
	
    public String getDescription() {return description;}
    public void setDescription(String str) {description = str;}
    
    /** The need for a short text to put on pull down menus made it necessary to add this function. By default
     * it takes the description field (possibly long) and will return a substring of it, indicating continuation with
     * with ellipsis (...)
     * @return Short description of this field
     */ 
    public String getShortDescription() { return short_description; }
    public void setShortDescription( String str ) { short_description = str; }
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