package etomica.action;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;

/**
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class BoxInflateAnisotropic extends BoxInflate{
    
    public BoxInflateAnisotropic(ISpace space){
    	super(space);
    	
        cVector = space.makeVector();
    	edgeVectorOld = new IVectorMutable[3];
        edgeVectorNew = new IVectorMutable[3];
        
        for(int i=0; i<edgeVectorOld.length; i++){
         	edgeVectorOld[i] = space.makeVector();
            edgeVectorNew[i] = space.makeVector();
        }
    }
    
    public BoxInflateAnisotropic(IBox box, ISpace space){
    	this(space);
    	setBox(box);
    }
    
    
    /**
     * Performs anisotropic inflation.
     */
    public void actionPerformed() {
        if(box == null) return;
        
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
    
        IVectorMutable translationVector = translator.getTranslationVector();
        
        for(int i=0; i<edgeVectorOld.length; i++){
         	edgeVectorOld[i].E(box.getBoundary().getEdgeVector(i));
        }

        double cx = cVector.getX(0);
        double deltacx = (cx-edgeVectorOld[2].getX(0));
        double cz = cVector.getX(2);
        double slope = deltacx/cz;
        
        //System.out.println("cVector: " + cVector.toString());
        
        IMoleculeList molecules = box.getMoleculeList();
        for(int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            translationVector.E(moleculeCenter.position(molecule));
            double h = translationVector.getX(2);
//            if(h<0.0){
//            	slope *= -1;
//            }
            scaleVector.E(new double[]{0.5*h*slope, 0.0, 0.0});
            System.out.println(" vector: " + translationVector.toString());
            System.out.println(h+ " scaleVector: " + scaleVector.toString());
            scaleVector.PE(-1.0);
            
            translationVector.TE(scaleVector);
            groupScaler.actionPerformed(molecule);
        }
//        System.exit(1);
        scaleVector.PE(1.0);

        // set the edgeVectors according to the scaling before passing it to BoundaryDeformablePeriodic
        // only scale the x-, y- and z-axes
        // for the x-component of edgeVector[2], it is not being scale. This is done so as to fluctuate
        // the beta-angle
        
        edgeVectorNew[0].E(edgeVectorOld[0]);
        edgeVectorNew[1].E(edgeVectorOld[1]);
        edgeVectorNew[2].E(cVector);
        
//        dimVector.E(box.getBoundary().getBoxSize());
//        System.out.println("dimVector: " + dimVector.toString());
//        dimVector.TE(scaleVector);
//        System.out.println("dimVector: " + dimVector.toString());
        
        ((BoundaryDeformablePeriodic)box.getBoundary()).setBoxSizeAngleFluctuation(edgeVectorNew);
    }
    
    
    public void undo(){
    	System.out.println("*******************************UNDO!!!!!!!!");
    	setCVector(edgeVectorOld[2]);
    	actionPerformed();
    }
        
    public void setCVector(IVectorMutable cVec){
    	cVector = cVec;
    }
    
    protected IVectorMutable[] edgeVectorOld, edgeVectorNew;
    protected IVectorMutable cVector;
	private static final long serialVersionUID = 1L;
}
