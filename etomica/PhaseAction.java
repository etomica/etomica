package etomica;

import java.awt.event.ActionEvent;

 /**
  * Superclass of classes that apply some elementary action (transformation) to a phase.
  * Inherits from javax.swing.AbstractAction, so may be registered as a listener to GUI components
  * 
  */
public abstract class PhaseAction extends etomica.Action {

    public static String getVersion() {return "PhaseAction:01.03.28.0/"+Action.VERSION;}

    protected Phase phase;
    public PhaseAction() {this(null);}
    public PhaseAction(Phase p) {
        super();
        phase = p;
    }
        
    public void setPhase(Phase p) {phase = p;}
    public Phase getPhase() {return phase;}
        
    public void actionPerformed(PhaseEvent pe) {
        actionPerformed(pe.getPhase());
    }
    
    public void actionPerformed() {actionPerformed(phase);}
    
    public abstract void actionPerformed(Phase p);

        
    //************* End of fields and methods for PhaseAction ************//
    
    //************* The remainder defines some subclasses of PhaseAction **********//
    
    public static final class Translate extends PhaseAction implements Action.Undoable {
        
        private Space.Vector translationVector;
        
        public Translate(Phase p) {
            super(p);
        }
        
        public void setTranslationVector(Space.Vector v) {translationVector.E(v);}
        public Space.Vector getTranslationVector() {return translationVector;}
        
        public void actionPerformed(Phase p) {doAction(p, translationVector);}
        
        public void actionPerformed(Space.Vector v) {doAction(phase, v);}
        
        public static void doAction(Phase p, Space.Vector v) {
            if(v == null || p == null) return;
            for(Molecule m=p.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.translateBy(v);
            }
        }
        public void attempt() {doAction(phase, translationVector);}
        public void undo() {
            translationVector.TE(-1);
            doAction(phase, translationVector);
            translationVector.TE(-1);
        }
    }//end of Translate
    
    public static final class ImposePbc extends PhaseAction implements Integrator.IntervalListener.ImposePbc {
        
        public ImposePbc(Phase p) {
            super(p);
        }
        
        public void actionPerformed(Phase p) {doAction(p);}
        
        public static void doAction(Phase p) {
            Space.Boundary boundary = p.boundary();
            for(Atom a = p.firstAtom(); a!=null; a=a.nextAtom()) {
                boundary.centralImage(a.coordinate);
            }
        }
        
        public void intervalAction(Integrator.IntervalEvent evt) {
            doAction(phase);
        }
        
        
    }//end of ImposePBC
    
    public static final class Inflate extends PhaseAction implements Action.Undoable {
        
        private double scale = 1.0;
        private transient Space.Vector temp;
        
        public Inflate(Phase p) {
            super(p);
            if(p != null) temp = p.parentSimulation().space().makeVector();
        }
                
        public void setScale(double s) {scale = s;}
        public double getScale() {return scale;}

        public void setPhase(Phase p) {
            super.setPhase(p);
            if(p != null && temp == null) temp = p.parentSimulation().space().makeVector();
        }
 
        public void actionPerformed(double s) {doAction(phase, s, temp);}
        public void actionPerformed(Phase p) {doAction(p, scale, temp);}

        //specific to 2D
        public static void doAction(Phase phase, double scale, int i, Space.Vector work){
            phase.boundary().dimensions().setComponent(i,scale*phase.boundary().dimensions().component(i));         
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                if (i==0){
                    work.setComponent(0, (scale-1.0)*((Space2D.Vector)m.position()).x);
                    work.setComponent(1,0.);
                }
                else {
                    work.setComponent(1, (scale-1.0)* ((Space2D.Vector)m.position()).y);
                    work.setComponent(0, 0.);                    
                }
                m.translateBy(work);   //displaceBy doesn't use temp
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
 //               if(display != null && i % 10 ==0) display.repaint();
            }
        }
        
        
        public static void doAction(Phase phase, double scale, Space.Vector work) {
            phase.boundary().inflate(scale);
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                work.Ea1Tv1(scale-1.0,m.position());
                m.translateBy(work); 
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
        
        public void attempt() {
            phase.boundary().inflate(scale);
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                temp.Ea1Tv1(scale-1.0,m.position());
                m.displaceBy(temp);   //displaceBy doesn't use temp
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
        public void undo() {
            phase.boundary().inflate(1.0/scale);
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.replace();
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
         
        public void attempt(int i){
            phase.boundary().dimensions().setComponent(i,scale*phase.boundary().dimensions().component(i));         
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                if (i==0){
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
 //               if(display != null && i % 10 ==0) display.repaint();
            }
        }
        public void undo(int i){
            phase.boundary().dimensions().setComponent(i,phase.boundary().dimensions().component(i)/scale);         
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.replace();
                m.atomIterator.reset();
                while(m.atomIterator.hasNext()) {phase.iteratorFactory().moveNotify(m.atomIterator.next());}
            }
        }
    }//end of Inflate 
}