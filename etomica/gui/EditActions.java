/**
 * EditActions
 *
 * The EditActions class is responsible for creating static action listeners to the Edit drop-down
 * menu JMenuItems of the EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package etomica.gui;

import etomica.beans.PhaseCustomizer;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.beans.*;
import javax.swing.JInternalFrame;

public class EditActions {
    /**
     * Static action listener that cuts a component from the FormDesign
     */
    public static final ActionListener CUT = new CutAction();

    /**
     * Static action listener that copies a component from the FormDesign
     */
    public static final ActionListener COPY = new CopyAction();
    
    /**
     * Static action listener that pastes a copied component to the FormDesign
     */
    public static final ActionListener PASTE = new PasteAction();
    
    /**
     * Static action listener that pastes a copied component to the FormDesign
     */
    public static final ActionListener CUSTOMIZE = new CustomizeAction();
    
    /**
     * Static action listener that displays information about the Etomica environment
     */
    public static final ActionListener REPORT= new ReportAction();
    
    /**
     * Static action listener that allows for two components to be bound to one another
     */
    public static final ActionListener BINDPROPERTY = new BindPropertyAction();
    
    /**
     * Static action listener that displays the current setup preferences of the Etomica environment
     */
    public static final ActionListener PREFERENCES = new PreferencesAction();

    private static Object obj = null;
    public static final void setObject(Object o) { obj = o; }
    public static final Object getObject() { return obj; }
    
    /**
     * This will eventually call a cut method for deleting an instance of a simulation component
     */
    private static class CutAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
                /*	    Wrapper wrapper = BeanBoxFrame.getCurrentWrapper();
	        Object bean = wrapper.getBean();
	        if (bean != BeanBoxFrame.getTopBox()) {
	            if (copy()) {
		            // succeeded in serializing the component
		            Container parent = wrapper.getParent();
		            BeanBoxFrame.setCurrentComponent(null);
		            if (parent != null) {
		                parent.remove(wrapper);
		            }
		            wrapper.cleanup();
	            }

	            bcss.remove(bean);
	        }
    	    
    */
        }// end of actionPerformed
    }// end of CutAction class
    
    /**
     * This will eventually call a copy method for making a copy of a simulation component
     */
    private static class CopyAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
                /*	    Wrapper wrapper = BeanBoxFrame.getCurrentWrapper();

	        BeanBoxFrame.setClipLabel(null);
	        BeanBoxFrame.setClipName(null);
	        try {
	            File f = new File(BeanBoxFrame.getClipFileName());
	            File dir = new File(f.getParent());
	            dir.mkdirs();
                FileOutputStream fos = new FileOutputStream(f);
                ObjectOutputStream oos = new ObjectOutputStream(fos);
	            // Ask the Wrapper to serialize the "naked" bean.
	            // This causes the Wrapper to remove all listeners
	            // before serializing, and add them back after.
	            wrapper.writeNakedBean(oos);
	            oos.close();
	            fos.close();
	            BeanBoxFrame.setClipLabel(wrapper.getBeanLabel());
	            BeanBoxFrame.setClipName(wrapper.getBeanName());	    
	            // We need to preserve information about whether or not a bean
	            // originated from a .ser file across copies and pastes.
	            // See paste(), below for more explanation.
	            BeanBoxFrame.setClipFromPrototypeInfo(wrapper.isFromPrototype());

	            pasteMenuItem.setEnabled(true);
	            return true;
	        } catch (Exception ex) {
	            error("Copy failed", ex);
	            pasteMenuItem.setEnabled(false);
	            return false;
	        }
    	    
    */
        }// end of actionPerformed
    }//end of CopyAction class
    
    /**
     * This will eventually call a paste method for pasting a copy of a simulation component to the FormDesign
     */
    private static class PasteAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
                /*	    synchronized (this) {
	            mouseClickEvent = null;
                }
	        try {
	            // Set the insert cursor before reading the clipboard.
	            setCursor(crosshairCursor);

	            SimpleClassLoader loader = SimpleClassLoader.ourLoader;
	            MyProducer p = new MyProducer(BeanBoxFrame.getClipFileName());
	            String clipName = BeanBoxFrame.getClipName();
	            String beanLabel = BeanBoxFrame.getClipLabel();

	            // We need to preserve information about whether or not a bean
	            // originated from a .ser file across copies and pastes.
        	    
	            // Beans that originate from a .ser file must be treated
	            // effectively as having hidden-state. As such, they must always
	            // be serialized. We need to propagate this information to 
	            // a newly instantiated beans' wrapper, when we do an insert
	            // into the beanbox, so that this information will be available
	            // later at code generation time.
	            boolean fromPrototypeInfo =BeanBoxFrame.getClipFromPrototypeInfo();
        		
	            Object bean = loader.instantiate(clipName, p);
	            doInsert(bean, beanLabel, clipName, true, fromPrototypeInfo);
	        } catch (Exception ex) {
	            error("Paste failed", ex);
	            pasteMenuItem.setEnabled(false);
	            setCursor(defaultCursor);
   	        }
        	    
    */
        }//end of actionPerformed
    }// end of PasteAction class

    /**
     * This will eventually call the Customizer for the selected bean
     */
    private static class CustomizeAction implements ActionListener {
        public void actionPerformed(ActionEvent event) {
            BeanInfo bi = null;
	        
	        try {
	            bi = Introspector.getBeanInfo(obj.getClass());
	        } 
	        catch (IntrospectionException ex) {
	            System.out.println("Couldn't introspect");
	        }
	        try {
	            Customizer c = ((Customizer)bi.getBeanDescriptor().getCustomizerClass().newInstance());
	            c.setObject(obj);
	            JInternalFrame frame = new javax.swing.JInternalFrame("Phase Customizer",true,true,true,true);
	            frame.setLocation(515, 60);
	            frame.setSize(((PhaseCustomizer)c).getPreferredSize());
	            frame.getContentPane().add((java.awt.Component)c);
	            Etomica.DesktopFrame.desktop.add(frame);
	            frame.setSelected(true);
            }
            catch (PropertyVetoException pve){}
            catch (InstantiationException ie){}
            catch (IllegalAccessException iae){}
        }//end of actionPerformed
    }// end of CustomizeAction class

    /**
     *  This will eventually call a reporting method for the current bean
     */
    private static class ReportAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
                /*	    Object bean = BeanBoxFrame.getCurrentBean();
	        if (bean == null) {
	            System.out.println("No current focus.");
	            return;
	        }
	        Report.report(bean.getClass());
    */
        }//end of actionPerformed
    }// end of ReportAction class

    /**
     * This will eventually call a bind property method for linking two simulation components
     */
    private static class BindPropertyAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }//end of actionPerformed
    }// end of BindPropertyAction class

    /**
     * Handles events from the Preferences MenuItem of the Edit Menu.  It displays the PreferenceFrame.
     */
    private static class PreferencesAction implements ActionListener {
                
        public void actionPerformed(ActionEvent event) {    
            new PreferenceFrame();
        }//end of actionPerformed
    }// end of PreferencesAction class
}// end of EditActions class
