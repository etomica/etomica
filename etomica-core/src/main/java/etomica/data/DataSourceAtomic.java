/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.atom.IAtom;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceAtomic {
    
    public IData getData(IAtom a);
    
    public IEtomicaDataInfo getAtomDataInfo();
    
    public DataTag getTag();
}
