package etomica;
import etomica.beans.*;
import etomica.units.Unit;
import java.awt.*;
import java.beans.*;

public class DisplayMeterHistogramBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Meter.class, MeterEditor.class);
        //PropertyEditorManager.registerEditor(Unit.class, DeviceUnitEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("label",      DisplayMeterHistogram.class),
                new PropertyDescriptor("makeResetButton",      DisplayMeterHistogram.class),
                new PropertyDescriptor("meter",      DisplayMeterHistogram.class),
                new PropertyDescriptor("name",       DisplayMeterHistogram.class),
                new PropertyDescriptor("updateInterval", DisplayMeterHistogram.class)};
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
}

