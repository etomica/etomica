package etomica.data;

import java.util.ArrayList;

import etomica.units.Dimension;

/**
 * Interface for a mutable DataInfo factory.  The factory is created based on 
 * an existing DataInfo object, but can be modified before being used to create
 * new DataInfo objects.
 */
public abstract class DataInfoFactory implements java.io.Serializable {
    /**
     * Creates a new instance using the info held by the template.  The 
     * template is not modified.
     */
    protected DataInfoFactory(DataInfo template) {
        label = template.getLabel();
        dimension = template.getDimension();
        tags = new ArrayList();
        Object[] templateTags = template.getTags();
        for (int i=0; i<templateTags.length; i++) {
            tags.add(templateTags[i]);
        }
    }
    
    /**
     * Creates a new DataInfo object using the information held by this factory.
     */
    public abstract DataInfo makeDataInfo();
    
    /**
     * Sets the label
     */
    public void setLabel(String newLabel) {
        label = newLabel;
    }
    
    /**
     * Returns the label
     */
    public String getLabel() {
        return label;
    }
    
    /**
     * Sets the dimension
     */
    public void setDimension(Dimension newDimension) {
        dimension = newDimension;
    }
    
    /**
     * Returns the dimension
     */
    public Dimension getDimension() {
        return dimension;
    }
    
    public ArrayList getTags() {
        return tags;
    }
    
    public void setTags(ArrayList newTags) {
        tags = (ArrayList)newTags.clone();
    }

    protected String label;
    protected Dimension dimension;
    protected ArrayList tags;
}