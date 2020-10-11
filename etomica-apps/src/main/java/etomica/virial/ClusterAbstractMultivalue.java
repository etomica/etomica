/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

public interface ClusterAbstractMultivalue extends ClusterAbstract {

    /**
     * @return number of values returned by getAllLastValues
     */
    int getNumValues();

    double[] getAllLastValues(BoxCluster box);
        
}
