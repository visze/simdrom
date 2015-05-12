package de.charite.compbio.simdrom.filter;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author Max Schubach <max.schubach@charite.de>
 */
public interface IFilter {
	
	public FilterType getFilterType();

	public VariantContext filter(VariantContext vc);

}