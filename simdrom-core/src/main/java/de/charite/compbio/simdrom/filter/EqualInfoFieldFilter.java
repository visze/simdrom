package de.charite.compbio.simdrom.filter;

public class EqualInfoFieldFilter extends AInfoFieldFilter {

	public EqualInfoFieldFilter(String info, Object type) {
		super(info, type);
	}

	@Override
	protected boolean compareInfoType(Object should, Object is) {
		return compare(should, is) == 0;
	}

}
