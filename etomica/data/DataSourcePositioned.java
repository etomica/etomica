/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.api.IBox;
import etomica.api.IVector;


/**
 * Interface for a DataSource that can return a value given an arbitrary
 * position.
 * 
 * @author Andrew Schultz
 */
public interface DataSourcePositioned {
    
    public void setBox(IBox box);
    
    public IData getData(IVector a);
    
    public IEtomicaDataInfo getPositionDataInfo();
    
    public DataTag getTag();
}