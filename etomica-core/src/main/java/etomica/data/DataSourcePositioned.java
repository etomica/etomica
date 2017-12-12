/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.box.Box;
import etomica.space.Vector;


/**
 * Interface for a DataSource that can return a value given an arbitrary
 * position.
 * 
 * @author Andrew Schultz
 */
public interface DataSourcePositioned {
    
    public void setBox(Box box);
    
    public IData getData(Vector a);
    
    public IDataInfo getPositionDataInfo();
    
    public DataTag getTag();
}
