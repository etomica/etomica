/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.data;

import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;

import java.util.ArrayList;
import java.util.List;

public class DataSourceAggregatorPotential extends DataSourceAggregator implements PotentialCallback {

    protected final PotentialCompute compute;
    protected final List<PotentialCallback> pairCallbacks;

    public DataSourceAggregatorPotential(PotentialCompute compute, IDataSourcePotential... dataSources) {
        super(dataSources);
        this.compute = compute;
        pairCallbacks = new ArrayList<>();
    }

    public IData getData() {
        boolean doForces = false;
        pairCallbacks.clear();
        for (IDataSource ds : dataSources) {
            IDataSourcePotential dsp = (IDataSourcePotential) ds;
            doForces = doForces || dsp.needsForces();
            if (dsp.needsPairCallback()) pairCallbacks.add(dsp.getPotentialCallback());
            dsp.doCallComputeAll(false);
        }
        compute.computeAll(doForces, pairCallbacks.size() > 0 ? this : null);
        return super.getData();
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        for (PotentialCallback pc : pairCallbacks) {
            pc.pairCompute(i, j, dr, u012);
        }
    }
}
