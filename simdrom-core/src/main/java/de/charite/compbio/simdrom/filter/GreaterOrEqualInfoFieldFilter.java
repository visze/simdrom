package de.charite.compbio.simdrom.filter;

public class GreaterOrEqualInfoFieldFilter extends AInfoFieldFilter {

	public GreaterOrEqualInfoFieldFilter(String info, Object type) {
		super(info, type);
	}

	@Override
	protected boolean compareInfoType(Object should, Object is) {
		return compare(should, is) == -1;
	}

}
