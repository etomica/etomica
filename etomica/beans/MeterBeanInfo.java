package etomica;
import etomica.*;
import java.awt.*;
import java.beans.*;
import java.lang.reflect.Method;

public class MeterBeanInfo extends SimpleBeanInfo {
    
    static {
//        PropertyEditorManager.registerEditor(Phase.class, PhaseEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            PropertyDescriptor phaseDescriptor =
                new PropertyDescriptor("phase", MeterAbstract.class);
//            phaseDescriptor.setPropertyEditorClass(PhaseEditor.class);

            return new PropertyDescriptor[] {
                phaseDescriptor,
                new PropertyDescriptor("mostRecent",Meter.class,"mostRecent",null),
                new PropertyDescriptor("average",Meter.class,"average",null),
                new PropertyDescriptor("error",Meter.class,"error",null),
                new PropertyDescriptor("active",MeterAbstract.class),
                new PropertyDescriptor("historying",MeterAbstract.class),
                new PropertyDescriptor("histogramming",MeterAbstract.class),
                new PropertyDescriptor("label",MeterAbstract.class),
                new PropertyDescriptor("updateInterval",MeterAbstract.class)
                };
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
    
//    public BeanDescriptor getBeanDescriptor() {
//        return new BeanDescriptor(Phase.class, PhaseCustomizer.class);
//    }
    
}

