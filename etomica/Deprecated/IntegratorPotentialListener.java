package simulate;
import java.util.*;

public interface IntegratorPotentialListener extends EventListener {
    public void integratorPotentialNotify(IntegratorPotentialEvent evt);
}