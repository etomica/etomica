package etomica.modifier;

import etomica.Modifier;
import etomica.units.Dimension;



/*
 * History
 * Created on Jan 31, 2005 by kofke
 */
class ModifierNMolecule extends Modifier {
        
        public ModifierNMolecule() {
            super(Dimension.QUANTITY);
        }
        
        public void setValue(double d) {
 //           if(initializing) return;
            if(d < 0) d = 0;
            boolean isPaused = this.selector.integrator.isPaused();
   //         if(!isPaused) {
                try {
                    this.selector.restartAction.actionPerformed();
                } catch (NullPointerException ex) {return;}
   //         }
            try {
                 this.selector.speciesAgent.setNMolecules((int)d);
            } catch(NullPointerException ex) {}
                try {
                    this.selector.restartAction.actionPerformed();
                } catch (NullPointerException ex) {return;}
            if(this.selector.display != null) this.selector.display.repaint();
            this.selector.integrator.reset();
        }
        public double getValue() {return (this.selector.speciesAgent!=null)?(double)this.selector.speciesAgent.moleculeCount():0;}
    }