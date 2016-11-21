package de.charite.compbio.simdrom.filter;

import java.util.Optional;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Interface of filters uised on the VCF files.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public interface IFilter {

	/**
	 * @return The {@link FilterType} of the filter.
	 */
	public FilterType getFilterType();

	/**
	 * Filter a variant by using this method. It can remove parts of the
	 * variant, or remove the variant completely. This returns <code>null</code>
	 * .
	 * 
	 * @param vc
	 *            variant to filter
	 * @return <code>null</code>, maybe party filtered or (if passed) the old
	 *         variant.
	 */
	public Optional<VariantContext> filter(VariantContext vc);
	
	public Optional<VariantContext> filter(Optional<VariantContext> vc);

}
