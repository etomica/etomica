/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.util.ArrayList;

import etomica.units.dimensions.Dimension;

/**
 * Interface for a mutable DataInfo factory.  The factory is created based on 
 * an existing DataInfo object, but can be modified before being used to create
 * new DataInfo objects.
 */
public abstract class DataInfoFactory implements java.io.Serializable, IEtomicaDataInfoFactory {

    /**
     * Creates a new instance using the info held by the template.  The 
     * template is not modified.
     */
    protected DataInfoFactory(IEtomicaDataInfo template) {
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
    public abstract IEtomicaDataInfo makeDataInfo();
    
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
    
    public DataTag[] getTags() {
        return (DataTag[])tags.toArray(new DataTag[tags.size()]);
    }
    
    public void setTags(DataTag[] newTags) {
        tags.clear();
        for (int i=0; i<newTags.length; i++) {
            tags.add(newTags[i]);
        }
    }

    protected String label;
    protected Dimension dimension;
    protected ArrayList tags;
}
