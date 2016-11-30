package de.charite.compbio.simdrom.filter;

public class LessOrEqualInfoFieldFilter extends AInfoFieldFilter {

	public LessOrEqualInfoFieldFilter(String info, Object type) {
		super(info, type);
	}

	@Override
	protected boolean compareInfoType(Object should, Object is) {
		int compareResult = compare(should, is);
		return compareResult == 1;
	}

}
