package etomica;
import etomica.*;
import java.awt.*;
import java.beans.*;

public class MeterDensityBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Phase.class, PhaseEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            PropertyDescriptor phaseDescriptor =
 //               new PropertyDescriptor("phase", 
 //                                      MeterDensity.class.getMethod("getPhase",null),
 //                                      MeterDensity.class.getMethod("setPhase",new Class[] {Phase.class}));
                new PropertyDescriptor("phase", 
                                       MeterDensity.class);
            phaseDescriptor.setPropertyEditorClass(PhaseEditor.class);
            return new PropertyDescriptor[] {
                phaseDescriptor,
                new PropertyDescriptor("active",MeterDensity.class),
                new PropertyDescriptor("historying",MeterDensity.class),
                new PropertyDescriptor("histogramming",MeterDensity.class),
                new PropertyDescriptor("label",MeterDensity.class),
                new PropertyDescriptor("updateInterval",MeterDensity.class)
                };
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
//        catch(NoSuchMethodException e) {
//            System.out.println("Error in MeterDensityBeanInfo: "+e);
//            return null;
//        }
    }
    
//    public BeanDescriptor getBeanDescriptor() {
//        return new BeanDescriptor(Phase.class, PhaseCustomizer.class);
//    }
    
    
    public Image getIcon(int iconType){
        String name = "";
        if(iconType == BeanInfo.ICON_COLOR_16x16){
            name="COLOR_16x16";
        }
        else if (iconType == BeanInfo.ICON_COLOR_32x32){
            name="COLOR_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_32x32){
            name="MONO_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_16x16){
            name="MONO_16x16";
        }
        else return null;
        return loadImage("Molecule_" + name + ".gif");
    }
}

