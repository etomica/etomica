/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.molecule.IMolecule;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceMolecular {
    
    public IData getData(IMolecule a);
    
    public IDataInfo getMoleculeDataInfo();
    
    public DataTag getTag();
}
