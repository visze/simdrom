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
	 * Filter a variant by using this method. It can remove parts of the variant, or remove the variant completely. Then
	 * {@link Optional#isPresent()} will be <code>false</code>.
	 * 
	 * @param vc
	 *            variant to filter
	 * @return An {@link Optional} including the variant. {@link Optional#isPresent()} will be <code>false</code> if the
	 *         filter removes the variant.
	 */
	public Optional<VariantContext> filter(VariantContext vc);

	/**
	 * Filter a variant (used ins an {@link Optional}) by using this method. It can remove parts of the variant, or
	 * remove the variant completely. Then {@link Optional#isPresent()} will be <code>false</code>.
	 * 
	 * @param vc
	 *            Variant included in an {@link Optional} to filter. If {@link Optional#isPresent()} <code>false</code>
	 *            nothing will be done.
	 * @return An {@link Optional} including the variant. {@link Optional#isPresent()} will be <code>false</code> if the
	 *         filter removes the variant.
	 */
	public Optional<VariantContext> filter(Optional<VariantContext> vc);

}
