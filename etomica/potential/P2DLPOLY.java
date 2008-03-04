package etomica.potential;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.WriteConfigurationP2DLPOLY;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.space.Space;

public class P2DLPOLY extends Potential{
	
	public P2DLPOLY(Space space){
		super(2, space);
       	configP2DLPOLY = new WriteConfigurationP2DLPOLY();
    	configP2DLPOLY.setConfName("CONFIG");

    
	}
	
	public double energy(IAtomSet atoms) {

		configP2DLPOLY.setMolecule(atoms.getAtom(0), atoms.getAtom(1));
		
		
		
		configP2DLPOLY.actionPerformed();
		int typeInteraction = configP2DLPOLY.getTypeInteraction();
		
		try{
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec("DLMULTI");
			int exitVal = proc.waitFor();
			FileReader fileReader = new FileReader("ConfigEnergy");			
			BufferedReader bufReader = new BufferedReader(fileReader);
			
			String line = bufReader.readLine();
			return Double.parseDouble(line);
			
			
		}catch (IOException e){
			throw new RuntimeException(e);
		}catch (InterruptedException err){
			throw new RuntimeException(err);
		}
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}


	public void setBox(IBox box) {
		configP2DLPOLY.setBox(box);
	}
	
	public WriteConfigurationP2DLPOLY getConfigP2DLPOLY() {
		return configP2DLPOLY;
	}

	public void setConfigP2DLPOLY(WriteConfigurationP2DLPOLY configP2DLPOLY) {
		this.configP2DLPOLY = configP2DLPOLY;
	}

	private IBox box;
	private WriteConfigurationP2DLPOLY configP2DLPOLY;

}
