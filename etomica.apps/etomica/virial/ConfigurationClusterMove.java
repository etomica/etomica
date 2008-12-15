package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class ConfigurationClusterMove extends ConfigurationCluster {

	public ConfigurationClusterMove(ISpace _space, IRandom random) {
		super(_space);
		this.random = random;
	}

	public void initializeCoordinates(IBox box) {
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		ClusterAbstract sampleCluster = clusterBox.getSampleCluster();
		while (sampleCluster.value(clusterBox)== 0){
		IAtomList list = box.getLeafList();
		for (int i=1;i<list.getAtomCount();i++){
			((IVectorRandom)((IAtomPositioned)list.getAtom(i)).getPosition()).setRandomInSphere(random);
			((IAtomPositioned)list.getAtom(i)).getPosition().TE(2);
			 clusterBox.trialNotify();
			 clusterBox.acceptNotify();
		}
		}
	}
   protected final IRandom random;
}
