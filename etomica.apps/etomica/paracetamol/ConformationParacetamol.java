package etomica.paracetamol;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/*
 * Geometry of Paracetamol Molecule
 * Constructed with angles and bond lengths
 * 
 * @author Tai Tan
 */

public class ConformationParacetamol extends Conformation{
	
	/*
	 * Bond Length [unit in Amstrom]
	 */
	
	// Benzene Ring 
	private double bondLengthC1C2 = 1.395;
	private double bondLengthC2C3 = 1.385;
	private double bondLengthC3C4 = 1.395;
	private double bondLengthC4C5 = 1.395;
	private double bondLengthC5C6 = 1.385;
	private double bondLengthC6C1 = 1.395; 
	
	private double bondLengthC4O1 = 1.352; //Benzene - OH
	private double bondLengthC1N1 = 1.394; //Benzene - N
	
	private double bondLengthN1C7 = 1.366; //N - C
	private double bondLengthC7O2 = 1.226; //C = O
	private double bondLengthC7C8 = 1.503; //C - C
	
	/*
	 * Bond Angles [unit is radian]
	 */
	
	private double angleC6C1C2 = 120.0*Math.PI/180;
	private double angleC2C3C4 = 120.0*Math.PI/180;
	private double angleC3C4C5 = 120.0*Math.PI/180;
	private double angleC5C6C1 = 120.0*Math.PI/180;
	
	private double angleO1C4C5 = 123.0*Math.PI/180; //Benzene - OH
	private double angleC6C1N1 = 118.0*Math.PI/180; //Benzene -N
	private double angleC1N1C7 = 129.0*Math.PI/180;	
	
	private final AtomIteratorArrayListSimple iterator;
	
	public ConformationParacetamol(Space space) {
		super(space);
		iterator = new AtomIteratorArrayListSimple();	
	}
	
	public void initializePositions(AtomSet List){
		
		iterator.setList(List);
		double x = 0.0;
		double y = 0.0;
		
		iterator.reset();
		
        IAtomPositioned c1 = (IAtomPositioned)iterator.nextAtom();
		c1.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c2 = (IAtomPositioned)iterator.nextAtom();
		x = x-bondLengthC1C2*Math.cos((angleC6C1C2)/2);
		y = y-bondLengthC1C2*Math.sin((angleC6C1C2)/2);
		c2.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c3 = (IAtomPositioned)iterator.nextAtom();
		x = x-bondLengthC2C3;
		c3.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c4 = (IAtomPositioned)iterator.nextAtom();
		x = x-bondLengthC3C4*Math.cos((angleC2C3C4/2));
		y = y+bondLengthC3C4*Math.sin((angleC2C3C4/2));
		c4.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned o1 = (IAtomPositioned)iterator.nextAtom();
		x = x-bondLengthC4O1*Math.cos((angleO1C4C5-120*Math.PI/180));
		y = y-bondLengthC4O1*Math.sin((angleO1C4C5-120*Math.PI/180));
		o1.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c5 = (IAtomPositioned)iterator.nextAtom();
		x = x+bondLengthC4O1*Math.cos((angleO1C4C5-120*Math.PI/180))+bondLengthC4C5*Math.cos((angleC3C4C5/2));
		y = y+bondLengthC4O1*Math.sin((angleO1C4C5-120*Math.PI/180))+bondLengthC4C5*Math.sin((angleC3C4C5/2));  
		c5.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c6 = (IAtomPositioned)iterator.nextAtom();
		x = x+bondLengthC5C6;
		c6.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned n1 = (IAtomPositioned)iterator.nextAtom();
		x = x+bondLengthC6C1*Math.cos((angleC5C6C1/2))+bondLengthC1N1*Math.cos((120*Math.PI/180-angleC6C1N1));
		y = y-bondLengthC6C1*Math.sin((angleC5C6C1/2))+bondLengthC1N1*Math.sin((120*Math.PI/180-angleC6C1N1));
		n1.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c7 = (IAtomPositioned)iterator.nextAtom();
		x = x+bondLengthN1C7*Math.cos((180*Math.PI/180-angleC1N1C7-2*Math.PI/180));
		y = y-bondLengthN1C7*Math.sin((180*Math.PI/180-angleC1N1C7-2*Math.PI/180));
		c7.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned c8 = (IAtomPositioned)iterator.nextAtom();
		x = x+bondLengthC7C8*Math.cos(17*Math.PI/180);
		y = y+bondLengthC7C8*Math.sin(17*Math.PI/180);
		c8.getPosition().E(new double [] {x, y, 0.0});
		
        IAtomPositioned o2 = (IAtomPositioned)iterator.nextAtom();
		x = x-bondLengthC7C8*Math.cos(17*Math.PI/180)-bondLengthC7O2*Math.cos(76*Math.PI/180);
		y = y-bondLengthC7C8*Math.sin(17*Math.PI/180)-bondLengthC7O2*Math.sin(76*Math.PI/180);
		o2.getPosition().E(new double [] {x, y, 0.0});
	
	}

	private static final long serialVersionUID = 1L;
}
