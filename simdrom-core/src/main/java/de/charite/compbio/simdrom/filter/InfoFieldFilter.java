package de.charite.compbio.simdrom.filter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

public class InfoFieldFilter implements IFilter {

	private final FilterType filterType = FilterType.INFO_FIELD_FILTER;

	private String info;
	private Object type;

	public InfoFieldFilter(String info, Object type) {
		this.info = info;
		this.type = type;
	}

	@Override
	public VariantContext filter(VariantContext vc) {
		CommonInfo infoField = vc.getCommonInfo();
		if (infoField.hasAttribute(info)) {
			Object val = infoField.getAttribute(info);
			if (val.getClass().isArray())
				return getVariantContextFromArray(vc, infoField.getAttribute(info));
			else if (val instanceof List)
				return getVariantContextFromArray(vc, ((List) val).toArray());
			else if (equalInfoType(val, type))
				return vc;
		}
		return null;
	}

	private VariantContext getVariantContextFromArray(VariantContext vc, Object infoField) {

		final int length = Array.getLength(infoField);
		if (length == vc.getAlternateAlleles().size()) {
			List<Allele> alleles = new ArrayList<Allele>();
			alleles.add(vc.getReference());
			for (int i = 0; i < length; i++) {
				if (equalInfoType(type, Array.get(infoField, i)))
					alleles.add(vc.getAlternateAllele(i));
			}
			if (alleles.size() <= 1)
				return null;
			else
				return new VariantContextBuilder(vc).alleles(alleles).make();
		} else {
			for (int i = 0; i < length; i++) {
				if (equalInfoType(type, Array.get(infoField, i)))
					return vc;
			}
		}
		return null;

	}

	private boolean equalInfoType(Object attribute, Object infoType) {
		return attribute.toString().equals(infoType.toString());
	}

	@Override
	public FilterType getFilterType() {
		return this.filterType;
	}
}
