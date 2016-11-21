package de.charite.compbio.simdrom.filter;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

/**
 * In implementation to use the VCF-Info column with a key {@link #info} and a value {@link #type} to filter out
 * variants that do not match to key=value.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class InfoFieldFilter extends AFilter {

	/**
	 * key in INFO column
	 */
	private String info;
	/**
	 * Value of the {@link #info} in INFO column
	 */
	private Object type;

	public InfoFieldFilter(String info, Object type) {
		super(FilterType.INFO_FIELD_FILTER);
		this.info = info;
		this.type = type;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.charite.compbio.simdrom.filter.IFilter#filter(htsjdk.variant. variantcontext.VariantContext)
	 */
	@Override
	public Optional<VariantContext> filter(Optional<VariantContext> optional_vc) {
		if (optional_vc.isPresent()) {
			VariantContext vc = optional_vc.get();
			CommonInfo infoField = vc.getCommonInfo();
			if (infoField.hasAttribute(info)) {
				Object val = infoField.getAttribute(info);
				if (val.getClass().isArray())
					return getVariantContextFromArray(vc, val);
				else if (val instanceof List)
					return getVariantContextFromArray(vc, ((List<?>) val).toArray());
				else if (equalInfoType(val, type))
					return optional_vc;
			}
		}
		return optional_vc;
	}

	private Optional<VariantContext> getVariantContextFromArray(VariantContext vc, Object infoField) {

		final int length = Array.getLength(infoField);

		// check if we have the same number of attributes than the same number
		// of ALT alleles.
		if (length == vc.getAlternateAlleles().size()) {
			List<Allele> alleles = new ArrayList<Allele>();
			alleles.add(vc.getReference());
			for (int i = 0; i < length; i++) {
				if (equalInfoType(type, Array.get(infoField, i)))
					alleles.add(vc.getAlternateAllele(i));
			}
			// no allele matches, return null
			if (alleles.size() <= 1)
				return Optional.empty();
			else
				return Optional.of(new VariantContextBuilder(vc).alleles(alleles).make());
		} else { // hack, no we add if one allele matches it.
			for (int i = 0; i < length; i++) {
				if (equalInfoType(type, Array.get(infoField, i)))
					return Optional.of(vc);
			}
		}
		return null;

	}

	private boolean equalInfoType(Object attribute, Object infoType) {
		return attribute.toString().equals(infoType.toString());
	}

}
