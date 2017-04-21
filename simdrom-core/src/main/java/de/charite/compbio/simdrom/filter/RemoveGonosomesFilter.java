/**
 * 
 */
package de.charite.compbio.simdrom.filter;

import java.util.Optional;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author max
 *
 */
public class RemoveGonosomesFilter extends AFilter {

	/**
	 * @param filterType
	 */
	public RemoveGonosomesFilter() {
		super(FilterType.GONOSOME_FILTER);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.charite.compbio.simdrom.filter.IFilter#filter(java.util.Optional)
	 */
	@Override
	public Optional<VariantContext> filter(Optional<VariantContext> vc) {
		if (vc.isPresent())
			if (vc.get().getContig()
					.matches("chr23|chrX|Chr23|ChrX|X|x|23|chrx|Chrx|chr24|chrY|Chr24|ChrY|Y|y|24|chry|Chry"))
				return Optional.empty();
		return vc;
	}

}
