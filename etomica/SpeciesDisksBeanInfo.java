package simulate;
import java.beans.*;

public class SpeciesDisksBeanInfo extends SimpleBeanInfo {

    public SpeciesDisksBeanInfo(){
        bd = new BeanDescriptor(SpeciesDisks.class, SpeciesDisksCustomizer.class);
      //  bd.setValue("hidden-state", Boolean.TRUE);
    }
    
    public BeanDescriptor getBeanDescriptor() {
        return bd;            
   }

   BeanDescriptor bd;
        }