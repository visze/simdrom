package de.charite.compbio.simdrom.filter;

import java.util.Optional;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Abstract filter class.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public abstract class AFilter implements IFilter {
	/**
	 * Filter type
	 */
	private final FilterType filterType;

	public AFilter(FilterType filterType) {
		this.filterType = filterType;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.charite.compbio.simdrom.filter.IFilter#filter(htsjdk.variant. variantcontext.VariantContext)
	 */
	@Override
	public Optional<VariantContext> filter(VariantContext vc) {
		Optional<VariantContext> optional_vc = Optional.of(vc);
		return filter(optional_vc);
	}

	@Override
	public FilterType getFilterType() {
		return this.filterType;
	}
}
