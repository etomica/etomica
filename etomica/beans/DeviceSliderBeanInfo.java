package etomica;
import etomica.beans.*;
import java.awt.*;
import java.beans.*;
import com.symantec.itools.vcafe.openapi.*;

public class DeviceSliderBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Modulator.class, ModulatorEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
//            PropertyDescriptor boundaryDescriptor =
//                new PropertyDescriptor("boundary", Phase.class);
//            boundaryDescriptor.setPropertyEditorClass(BoundaryEditor.class);
//            IndexedPropertyDescriptor moleculePositionDescriptor =
//                new IndexedPropertyDescriptor("moleculePosition", Phase.class);
//            moleculePositionDescriptor.setPropertyEditorClass(MoleculePositionEditor.class);
//            PropertyDescriptor meterListDescriptor =
//                new PropertyDescriptor("meterList", Phase.class);
//            meterListDescriptor.setPropertyEditorClass(MeterListEditor.class);
//            return new PropertyDescriptor[] {boundaryDescriptor, meterListDescriptor};
            PropertyDescriptor sliderDescriptor = new PropertyDescriptor("slider", DeviceSlider.class);
            sliderDescriptor.setConstrained(true);
            return new PropertyDescriptor[] {
                new PropertyDescriptor("modulator", DeviceSlider.class),
                new PropertyDescriptor("minimum", DeviceSlider.class),
                new PropertyDescriptor("maximum", DeviceSlider.class),
                sliderDescriptor
            };
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
}

