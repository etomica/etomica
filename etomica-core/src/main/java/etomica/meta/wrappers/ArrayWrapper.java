package etomica.meta.wrappers;

import etomica.meta.WrapperIndex;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ArrayWrapper extends CollectionWrapper<Object[]> {

    public ArrayWrapper(Object[] wrapped) {
        super(wrapped);
    }

    @Override
    public List<Wrapper<?>> getElements() {
        List<Wrapper<?>> elements = new ArrayList<>();
        for(Object el : wrapped) {
            Wrapper<?> wrapper = WrapperIndex.getWrapper(el);
            elements.add(wrapper);
        }

        return elements;
    }
}
