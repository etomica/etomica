/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleBDArray;
import etomica.data.types.DataDoubleBDArray.DataInfoDoubleBDArray;
import etomica.units.dimensions.Null;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * This meter handles direct sampling with BigDecimal precision for averages.
 */
public class MeterVirialBD implements IDataSource {

    protected static final BigDecimal BDZERO = new BigDecimal(0);

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialBD(ClusterAbstract[] aClusters) {
		clusters = aClusters;
		mc = new MathContext(40);
        data = new DataDoubleBDArray(1, 40);
        dataInfo = new DataInfoDoubleBDArray("Cluster Value",Null.DIMENSION, new int[]{1}, 40);
        tag = new DataTag();
        dataInfo.addTag(tag);
        DataDoubleBDArray.addBDZero(BDZERO);
	}

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
        
        BigDecimal[] x = data.getData();
        double v = clusters[1].value(box);
        if (v==0) {
            x[0] = BDZERO;
            return data;
        }
        double pi = box.getSampleCluster().value(box);

        if (pi == 0 || pi == Double.POSITIVE_INFINITY || Double.isNaN(pi)) throw new RuntimeException("oops "+pi);
        
        x[0] = new BigDecimal(v, mc).divide(new BigDecimal(pi, mc), mc);
        return data;
    }
    
    public ClusterAbstract[] getClusters() {
        return clusters;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }
    
    protected final ClusterAbstract clusters[];
	private final DataDoubleBDArray data;
	private final MathContext mc;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
}
