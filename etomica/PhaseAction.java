package etomica;

import java.awt.event.ActionEvent;

 /**
  * Superclass of classes that apply some elementary action (transformation) to a phase.
  * Inherits from javax.swing.AbstractAction, so may be registered as a listener to GUI components
  * 
  */
public abstract class PhaseAction extends etomica.Action {

    protected Phase phase;
    private DisplayPhase display; 
    public PhaseAction() {this(null);}
    public PhaseAction(Phase p) {
        super();
        phase = p;
    }
        
    public void setPhase(Phase p) {phase = p;}
    public Phase getPhase() {return phase;}
    
    
    public void actionPerformed(Phase p) {
        setPhase(p);
        actionPerformed();
    }
    
    public abstract void actionPerformed();
    
        
    //************* End of fields and methods for PhaseAction ************//
    
    //************* The remainder defines some subclasses of PhaseAction **********//
    
    public static final class Inflate extends PhaseAction implements Action.Retractable {
        
        private double scale = 1.0;
        private transient Space.Vector temp;
        
          public Inflate() {
          this(null);
        }
        public Inflate(Phase p) {
            super(p);
            temp = p.parentSimulation().space().makeVector();
        }
 
        public void actionPerformed(double s) {
            setScale(s);
            actionPerformed();
        }
        public void actionPerformed(Phase p, double s) {
            setScale(s);
            this.setPhase(p);
            actionPerformed();
        }
        public void actionPerformed(Phase p,  double s, int i){
          
           setScale(s);
           this.setPhase(p); 
             
            phase.boundary().dimensions().setComponent(i,s*phase.boundary().dimensions().component(i));         
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                if ( i==0){
                    temp.setComponent(0, (scale-1.0)*((Space2D.Vector)m.position()).x);
                    temp.setComponent(1,0.);
                }
                else {
                    temp.setComponent(1, (scale-1.0)* ((Space2D.Vector)m.position()).y);
                    temp.setComponent(0, 0.);                    
                }
                m.displaceBy(temp);   //displaceBy doesn't use temp
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
                if(display != null && i % 10 ==0) display.repaint();
            }
        }
        public void actionPerformed() {
            phase.boundary().inflate(scale);
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                temp.Ea1Tv1(scale-1.0,m.position());
                m.displaceBy(temp);   //displaceBy doesn't use temp
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
        public void retractAction() {
            phase.boundary().inflate(1.0/scale);
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.replace();
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
         
        public void retractAction(int i){
            phase.boundary().dimensions().setComponent(i,phase.boundary().dimensions().component(i)/scale);         
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.replace();
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
        public void setScale(double s) {scale = s;}
        public double getScale() {return scale;}
    }//end of Inflate 
}