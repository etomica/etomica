package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.space.ISpace;

public class ConfigurationClusterChain extends ConfigurationCluster {

	public ConfigurationClusterChain(ISpace _space) {
		super(_space);
	}

	public void initializeCoordinates(IBox box) {
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IAtomList list = box.getLeafList();
		for (int i=1;i<list.getAtomCount();i++){
			((IVectorMutable)list.getAtom(i).getPosition()).setX(0, 0.9*i);
		 }
		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 //System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
}
