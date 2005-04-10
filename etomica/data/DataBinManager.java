package etomica.data;

/**
 * Interfaces for class that manages one or more DataBin instances and initiates
 * some action when their data are changed.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public interface DataBinManager {

    /**
     * Method called by DataBin when its data has been changed.
     * 
     * @param sourceBin
     *            the bin that is notifying of the change
     */
    public void dataChangeNotify(DataBin sourceBin);

}